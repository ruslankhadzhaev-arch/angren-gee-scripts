// ==========================================================
//    ANGREN: PHD BENCHMARK EDITION (1990-2024)
//    Methodology:
//    1. Projection‑Consistent Area Coverage (RefProj Locked)
//    2. Optimized Two-Step Filtering
//    3. Full Metadata & Null-Safety
// ==========================================================

// --- НАСТРОЙКИ ---
var MIN_PX  = 10;  // Порог валидности сцены (>= 10 px)
var S_START = 6;
var S_END   = 8;  // Август

// 1. ROI (~6×6 км)
var roi = ee.Geometry.Rectangle([70.1446, 41.0091, 70.2158, 41.0632]);
Map.centerObject(roi, 13);
Map.addLayer(roi, {color: 'red'}, 'ROI Angren');

// ==========================================================
// 2. ФУНКЦИИ
// ==========================================================
function getCloudMask(img) {
  var qa = img.select('QA_PIXEL');
  return qa.bitwiseAnd(1 << 0).eq(0)
    .and(qa.bitwiseAnd(1 << 1).eq(0))
    .and(qa.bitwiseAnd(1 << 2).eq(0))
    .and(qa.bitwiseAnd(1 << 3).eq(0))
    .and(qa.bitwiseAnd(1 << 4).eq(0))
    .and(qa.bitwiseAnd(1 << 5).eq(0));
}

function processTemp(img, bandName) {
  return img.select(bandName)
           .multiply(0.00341802).add(149.0).subtract(273.15)
           .rename('LST');
}

function calcIndices(img) {
  var ndvi = img.normalizedDifference(['NIR', 'Red']).rename('NDVI');
  var savi = img.expression(
    '((NIR - RED) / (NIR + RED + 0.5)) * 1.5',
    {'NIR': img.select('NIR'), 'RED': img.select('Red')}
  ).rename('SAVI');
  return img.addBands([ndvi, savi]);
}

// ==========================================================
// 3. СБОР LANDSAT-КОЛЛЕКЦИЙ
// ==========================================================
var oldSensors = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
  .merge(ee.ImageCollection('LANDSAT/LE07/C02/T1_L2'))
  .filterBounds(roi)
  .filter(ee.Filter.calendarRange(S_START, S_END, 'month'))
  .map(function(img) {
    var mask = getCloudMask(img);
    var optical = img.select(['SR_B3', 'SR_B4']).updateMask(mask)
       .multiply(0.0000275).add(-0.2).rename(['Red','NIR']);
    var temp    = processTemp(img, 'ST_B6').updateMask(mask);
    return optical.addBands(temp).set('system:time_start', img.get('system:time_start'));
  });

var newSensors = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2'))
  .filterBounds(roi)
  .filter(ee.Filter.calendarRange(S_START, S_END, 'month'))
  .map(function(img) {
    var mask = getCloudMask(img);
    var optical = img.select(['SR_B4', 'SR_B5']).updateMask(mask)
       .multiply(0.0000275).add(-0.2).rename(['Red','NIR']);
    var temp    = processTemp(img, 'ST_B10').updateMask(mask);
    return optical.addBands(temp).set('system:time_start', img.get('system:time_start'));
  });

var landsatFull = oldSensors.merge(newSensors).sort('system:time_start');

// ERA5 (Климат)
var era5 = ee.ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR')
  .filterBounds(roi)
  .filter(ee.Filter.calendarRange(S_START, S_END, 'month'))
  .select(['total_precipitation_sum','volumetric_soil_water_layer_1']);

// ==========================================================
// 4. ЭТАЛОННАЯ ПРОЕКЦИЯ И ПЛОЩАДИ
// ==========================================================
var refProj  = ee.Image(landsatFull.first()).select('Red').projection();
var roiAreaGeo = ee.Number(roi.area(1));

// Знаменатель – сумма площадей пикселей 30 м в той же проекции
var roiArea30m = ee.Number(
  ee.Image.pixelArea()
    .reduceRegion({
      reducer: ee.Reducer.sum(),
      geometry: roi,
      scale: 30,
      crs: refProj,        // фиксируем CRS
      maxPixels: 1e9,
      tileScale: 4
    }).get('area')
);

// ==========================================================
// 5. ГЛАВНЫЙ ЦИКЛ (1990–2024)
// ==========================================================
var years = ee.List.sequence(1990, 2024);

var finalData = ee.FeatureCollection(years.map(function(y) {
  // A. Климат
  var eraYear = era5.filter(ee.Filter.calendarRange(y, y, 'year'));
  var precip = eraYear.select('total_precipitation_sum').sum()
      .reduceRegion({reducer: ee.Reducer.mean(), geometry: roi, scale: 11000})
      .get('total_precipitation_sum');
  precip = ee.Number(ee.Algorithms.If(precip, precip, 0));
  var precip_mm = precip.multiply(1000);

  var ssm = eraYear.select('volumetric_soil_water_layer_1').mean()
      .reduceRegion({reducer: ee.Reducer.mean(), geometry: roi, scale: 11000})
      .get('volumetric_soil_water_layer_1');
  ssm = ee.Number(ee.Algorithms.If(ssm, ssm, 0));

  // B. Landsat: фильтрация по количеству пикселей (два шага)
  var landsatFiltered = landsatFull.filter(ee.Filter.calendarRange(y, y, 'year'));

  // Шаг 1: считаем valid_pixels
  var withPx = landsatFiltered.map(function(img) {
    var pxCount = img.select('Red').reduceRegion({
      reducer: ee.Reducer.count(),
      geometry: roi,
      scale: 30,
      maxPixels: 1e9,
      tileScale: 4
    }).get('Red');
    pxCount = ee.Number(ee.Algorithms.If(pxCount, pxCount, 0));
    return img.set('valid_pixels', pxCount);
  });

  // фильтр по MIN_PX
  var validImages = withPx.filter(ee.Filter.gte('valid_pixels', MIN_PX));
  var count = validImages.size();

  // Шаг 2: считаем покрытие по площади для каждой сцены
  validImages = validImages.map(function(img) {
    var sceneArea = ee.Image.pixelArea()
      .updateMask(img.select('Red').mask())
      .reduceRegion({
        reducer: ee.Reducer.sum(),
        geometry: roi,
        scale: 30,
        crs: refProj,
        maxPixels: 1e9,
        tileScale: 4
      }).get('area');
    sceneArea = ee.Number(ee.Algorithms.If(sceneArea, sceneArea, 0));
    var ratio = sceneArea.divide(roiArea30m).min(1).max(0);
    return img.set('scene_cov_area_ratio', ratio);
  });

  // Средняя площадь покрытия сцен (0..1)
  var covSceneAreaMean = ee.Number(ee.Algorithms.If(
    count.gt(0),
    validImages.aggregate_mean('scene_cov_area_ratio'),
    0
  ));

  return ee.Algorithms.If(
    count.gt(0),
    // Есть валидные снимки
    (function() {
      // Медианный композит
      var composite = validImages.select(['Red','NIR','LST']).median();

      // Покрытие композита по площади
      var compValidArea = ee.Image.pixelArea()
        .updateMask(composite.select('Red').mask())
        .reduceRegion({
          reducer: ee.Reducer.sum(),
          geometry: roi,
          scale: 30,
          crs: refProj,
          maxPixels: 1e9,
          tileScale: 4
        }).get('area');
      compValidArea = ee.Number(ee.Algorithms.If(compValidArea, compValidArea, 0));
      var covCompAreaRatio = compValidArea.divide(roiArea30m).min(1).max(0);

      // Индексы и статистика
      var compositeIndices = calcIndices(composite);
      var reducer = ee.Reducer.mean().combine({
        reducer2: ee.Reducer.stdDev(),
        sharedInputs: true
      });
      var stats = compositeIndices.reduceRegion({
        reducer: reducer,
        geometry: roi,
        scale: 30,
        maxPixels: 1e9,
        tileScale: 4
      });

      return ee.Feature(null, {
        'Year': y,
        // Метаданные
        'Has_Landsat': ee.Number(1),
        'MIN_PX': MIN_PX,
        'SeasonStart': S_START,
        'SeasonEnd': S_END,
        'ROI_Area_Geo_m2': roiAreaGeo,
        'ROI_Area_30m_m2': roiArea30m,
        'Coverage_Method': 'pixelArea(RefProj)/pixelArea(RefProj)',
        // Основные данные
        'NDVI_mean': stats.get('NDVI_mean'), 'NDVI_std': stats.get('NDVI_stdDev'),
        'SAVI_mean': stats.get('SAVI_mean'), 'SAVI_std': stats.get('SAVI_stdDev'),
        'LST_mean_C': stats.get('LST_mean'), 'LST_std_C': stats.get('LST_stdDev'),
        'PP_mm': precip_mm,
        'SSM_m3m3': ssm,
        // Качество
        'Image_Count': count,
        'Cov_Scene_Mean': covSceneAreaMean,
        'Cov_Composite': covCompAreaRatio
      });
    })(),
    // Нет валидных снимков
    ee.Feature(null, {
      'Year': y,
      'Has_Landsat': ee.Number(0),
      'MIN_PX': MIN_PX,
      'SeasonStart': S_START,
      'SeasonEnd': S_END,
      'ROI_Area_Geo_m2': roiAreaGeo,
      'ROI_Area_30m_m2': roiArea30m,
      'Coverage_Method': 'pixelArea(RefProj)/pixelArea(RefProj)',
      'NDVI_mean': null, 'NDVI_std': null,
      'SAVI_mean': null, 'SAVI_std': null,
      'LST_mean_C': null, 'LST_std_C': null,
      'PP_mm': precip_mm,
      'SSM_m3m3': ssm,
      'Image_Count': 0,
      'Cov_Scene_Mean': 0,
      'Cov_Composite': 0
    })
  );
}));

finalData = finalData.sort('Year');

// ==========================================================
// 6. ГРАФИЧЕСКИЙ АНАЛИЗ
// ==========================================================

print('BENCHMARK ANALYSIS COMPLETE.');

print('Critical gaps (Composite < 0.3):', finalData.filter(ee.Filter.lt('Cov_Composite', 0.3)));

var chartClimate = ui.Chart.feature.byFeature(finalData, 'Year', ['PP_mm','SSM_m3m3'])
  .setChartType('ComboChart')
  .setOptions({
    title: 'Climate Context (ERA5)',
    seriesType: 'line',
    series: {0:{type:'bars', targetAxisIndex:0, color:'#1f77b4'},
             1:{type:'line', targetAxisIndex:1, color:'#8c564b', lineWidth:3}},
    vAxes: {0: {title:'Precip (mm)'}, 1:{title:'SSM (m³/m³)'}}
  });
print(chartClimate);

var chartVeg = ui.Chart.feature.byFeature(finalData, 'Year', ['SAVI_mean','SSM_m3m3'])
  .setChartType('ComboChart')
  .setOptions({
    title: 'Vegetation Response: SAVI vs Soil Moisture',
    seriesType: 'line',
    series: {0:{targetAxisIndex:0, color:'green', lineWidth:3},
             1:{targetAxisIndex:1, color:'brown', lineDashStyle:[4,4]}},
    vAxes: {0:{title:'SAVI'}, 1:{title:'SSM (m³/m³)'}}
  });
print(chartVeg);

var chartLST = ui.Chart.feature.byFeature(finalData, 'Year', ['LST_mean_C'])
  .setChartType('LineChart')
  .setOptions({
    title: 'Surface Temperature Trend (LST)',
    colors: ['red'],
    lineWidth:3,
    pointSize:4,
    vAxis: {title:'Temperature (°C)'},
    hAxis: {title:'Year', format:'yyyy'},
    trendlines: {0:{color:'black', opacity:0.5, showR2:true, visibleInLegend:true}}
  });
print(chartLST);

var chartQual = ui.Chart.feature.byFeature(finalData, 'Year', ['Cov_Scene_Mean','Cov_Composite'])
  .setChartType('LineChart')
  .setOptions({
    title: 'Data Reliability: Scene vs Composite Area Coverage',
    legend: {position:'bottom'},
    series: {0:{color:'#aaaaaa', lineDashStyle:[2,2], lineWidth:2},
             1:{color:'#000000', lineWidth:3}},
    vAxis: {title:'Coverage Area Ratio (0–1)', viewWindow:{min:0, max:1.1}},
    hAxis: {title:'Year', format:'yyyy'}
  });
print(chartQual);

var chartCount = ui.Chart.feature.byFeature(finalData, 'Year', ['Image_Count'])
  .setChartType('ColumnChart')
  .setOptions({
    title: 'Data Availability: Valid Landsat Images',
    legend:{position:'none'},
    colors:['#4682B4'],
    vAxis:{title:'Count'},
    hAxis:{title:'Year', format:'yyyy'}
  });
print(chartCount);

// ==========================================================
// 7. ЭКСПОРТ
// ==========================================================
Export.table.toDrive({
  collection: finalData,
  description: 'Angren_PhD_Benchmark_Dataset',
  fileFormat: 'CSV',
  selectors: [
    'Year','Has_Landsat',
    'MIN_PX','SeasonStart','SeasonEnd','Coverage_Method',
    'ROI_Area_Geo_m2','ROI_Area_30m_m2',
    'NDVI_mean','NDVI_std',
    'SAVI_mean','SAVI_std',
    'LST_mean_C','LST_std_C',
    'PP_mm','SSM_m3m3',
    'Image_Count','Cov_Scene_Mean','Cov_Composite'
  ]
});
