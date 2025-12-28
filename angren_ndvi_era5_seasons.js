/**************************************
 * Angren NDVI + ERA5 — сезоны 2019–2024
 * Landsat + Sentinel-2 + ERA5-Land
 * Таблица: год, сезон, NDVI-статистика и климат
 **************************************/

// --------- 0. ROI (6×6 км, как в таблицах NDVI) ---------
var roi = ee.Geometry.Rectangle([
  70.1446, 41.0091,   // lon_min, lat_min
  70.2158, 41.0632    // lon_max, lat_max
]);

Map.centerObject(roi, 12);
Map.addLayer(roi, {color: 'red'}, 'ROI 6×6 km', false);

// Период анализа (совпадает с ERA5 и S2)
var startYear = 2019;
var endYear   = 2024;
var years = ee.List.sequence(startYear, endYear);

// Сезоны (как ты использовал для ERA5: весна, лето, осень)
var seasons = ee.List([
  {name: 'spring', start: 3, end: 5},
  {name: 'summer', start: 6, end: 8},
  {name: 'autumn', start: 9, end: 10}
]);

// ================== 1. МАСКИ И NDVI ДЛЯ LANDSAT И S2 ==================

/** Маска для Landsat L2 (QA_PIXEL): убирает тени, снег, облака, воду */
function maskLandsatL2(img) {
  var qa = img.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(1 << 3).eq(0)   // cloud shadow
    .and(qa.bitwiseAnd(1 << 4).eq(0))      // snow/ice
    .and(qa.bitwiseAnd(1 << 5).eq(0))      // cloud
    .and(qa.bitwiseAnd(1 << 7).eq(0));     // water
  return img.updateMask(mask);
}

/** Добавляет NDVI и NBR к сценам Landsat (LT05/LE07/LC08/LC09) */
function addIndicesLandsat(img) {
  var hasB5 = img.bandNames().contains('SR_B5'); // для LC08/LC09

  // NDVI
  var ndvi = ee.Image(ee.Algorithms.If(
    hasB5,
    img.normalizedDifference(['SR_B5','SR_B4']).rename('NDVI'),
    img.normalizedDifference(['SR_B4','SR_B3']).rename('NDVI')
  ));

  // NBR (нам сейчас не нужен, но оставим для совместимости)
  var nbr  = ee.Image(ee.Algorithms.If(
    hasB5,
    img.normalizedDifference(['SR_B5','SR_B7']).rename('NBR'),
    img.normalizedDifference(['SR_B4','SR_B7']).rename('NBR')
  ));

  return img.addBands([ndvi, nbr]);
}

/** Маска для Sentinel-2: растительность (4) + голая почва (5) */
function maskS2(img){
  var scl = img.select('SCL');
  var good = scl.eq(4).or(scl.eq(5)); // vegetation + bare soil
  var bad  = scl.eq(3).or(scl.eq(8)).or(scl.eq(9)).or(scl.eq(10)); // тени, облака, снег
  return img.updateMask(good.and(bad.not()));
}

/** Добавляет NDVI и NBR к сценам Sentinel-2 SR */
function addIndicesS2(img) {
  var ndvi = img.normalizedDifference(['B8','B4']).rename('NDVI');
  var nbr  = img.normalizedDifference(['B8','B12']).rename('NBR');
  return img.addBands([ndvi, nbr]);
}

// ---------- Коллекции NDVI ----------

// Landsat 8 + 9 L2 (2013–…); для 2019–2024 этого достаточно
var landsatAll = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2'))
  .filterBounds(roi)
  .filterDate(startYear + '-01-01', (endYear + 1) + '-01-01')
  .map(maskLandsatL2)
  .map(addIndicesLandsat);

// Sentinel-2 SR (можно заменить на S2_SR_HARMONIZED при желании)
var S2 = ee.ImageCollection('COPERNICUS/S2_SR')
  .filterBounds(roi)
  .filterDate(startYear + '-01-01', (endYear + 1) + '-01-01')
  .filter(ee.Filter.lte('CLOUDY_PIXEL_PERCENTAGE', 40))
  .map(maskS2)
  .map(addIndicesS2);

print('Landsat сцены (2019–2024):', landsatAll.size());
print('Sentinel-2 сцены (2019–2024):', S2.size());

// ================== 2. ERA5-LAND (осадки и влажность почвы) ==================

var era = ee.ImageCollection('ECMWF/ERA5_LAND/HOURLY')
  .filterBounds(roi)
  .filterDate(startYear + '-01-01', (endYear + 1) + '-01-01');

// НОВАЯ КОЛЛЕКЦИЯ СУТОЧНЫХ АГРЕГАТОВ
var eraDaily = ee.ImageCollection('ECMWF/ERA5_LAND/DAILY_AGGR')
  .filterBounds(roi)
  .filterDate(startYear + '-01-01', (endYear + 1) + '-01-01');

// НОВАЯ ФУНКЦИЯ ERA5 ПО СЕЗОНАМ
function getEraSeasonStats(year, startMonth, endMonth, seasonName) {
  year = ee.Number(year);
  var start = ee.Date.fromYMD(year, startMonth, 1);
  var end   = ee.Date.fromYMD(year, endMonth, 1).advance(1, 'month');

  var collDaily = eraDaily.filterDate(start, end);

  // суточные суммы уже посчитаны, просто суммируем дни
  var tpSum = collDaily.select('total_precipitation_sum')
      .sum()
      .multiply(1000)          // м → мм
      .rename('tp_mm');

  // влажность почвы оставляем из HOURLY
  var collHourly = era.filterDate(start, end);
  var swvl1Mean = collHourly.select('volumetric_soil_water_layer_1')
      .mean()
      .rename('swvl1');

  var img = tpSum.addBands(swvl1Mean);

  var stats = img.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: roi,
    scale: 9000,
    maxPixels: 1e9
  });

  return ee.Feature(null, {
    year  : year,
    season: seasonName,
    tp_mm : stats.get('tp_mm'),
    swvl1 : stats.get('swvl1')
  });
}


// ================== 3. СЕЗОННЫЕ СТАТИСТИКИ NDVI ==================

/** Универсальная функция: NDVI-коллекция -> сезонные статистики по ROI */
function getNdviSeasonStats(collection, year, startMonth, endMonth, seasonName, sensorName, scale) {
  year = ee.Number(year);

  var start = ee.Date.fromYMD(year, startMonth, 1);
  var end   = ee.Date.fromYMD(year, endMonth, 1).advance(1, 'month');

  var coll = collection
      .filterDate(start, end)
      .filterBounds(roi);

  // сезонный композит NDVI
  var ndviSeason = coll.select('NDVI').mean().rename('ndvi');

  var stats = ndviSeason.reduceRegion({
    reducer: ee.Reducer.minMax()
      .combine(ee.Reducer.mean(), '', true)
      .combine(ee.Reducer.stdDev(), '', true),
    geometry: roi,
    scale: scale,
    maxPixels: 1e9
  });

  return ee.Feature(null, {
    'year'        : year,
    'season'      : seasonName,
    'sensor'      : sensorName,
    'ndvi_min'    : stats.get('ndvi_min'),
    'ndvi_max'    : stats.get('ndvi_max'),
    'ndvi_mean'   : stats.get('ndvi_mean'),
    'ndvi_std'    : stats.get('ndvi_stdDev')
  });
}

// --------- Объединённая функция: Landsat + S2 + ERA5 для одного года/сезона ---------
function seasonCombinedStats(year, seasonObj) {
  year = ee.Number(year);
  var s = ee.Dictionary(seasonObj);
  var sName = ee.String(s.get('name'));
  var startM = ee.Number(s.get('start'));
  var endM   = ee.Number(s.get('end'));

  // NDVI по Landsat
  var lsFeat = getNdviSeasonStats(
    landsatAll, year, startM, endM, sName, 'Landsat', 30
  );

  // NDVI по Sentinel-2
  var s2Feat = getNdviSeasonStats(
    S2, year, startM, endM, sName, 'Sentinel-2', 10
  );

  // ERA5-Land
  var eraFeat = getEraSeasonStats(year, startM, endM, sName);

  return ee.Feature(null, {
    'year'        : year,
    'season'      : sName,

    // Landsat
    'ndvi_min_ls'  : lsFeat.get('ndvi_min'),
    'ndvi_max_ls'  : lsFeat.get('ndvi_max'),
    'ndvi_mean_ls' : lsFeat.get('ndvi_mean'),
    'ndvi_std_ls'  : lsFeat.get('ndvi_std'),

    // Sentinel-2
    'ndvi_min_s2'  : s2Feat.get('ndvi_min'),
    'ndvi_max_s2'  : s2Feat.get('ndvi_max'),
    'ndvi_mean_s2' : s2Feat.get('ndvi_mean'),
    'ndvi_std_s2'  : s2Feat.get('ndvi_std'),

    // ERA5-Land
    'tp_mm'        : eraFeat.get('tp_mm'),
    'swvl1'        : eraFeat.get('swvl1')
  });
}

// --------- Строим таблицу для всех лет и сезонов ---------
var features = years.map(function(y) {
  y = ee.Number(y);
  return seasons.map(function(s) {
    return seasonCombinedStats(y, s);
  });
}).flatten();

var finalStats = ee.FeatureCollection(features);
print('Сезонные NDVI (Landsat+S2) + ERA5 2019–2024', finalStats);

// ================== 4. ЭКСПОРТ ТАБЛИЦЫ В EXCEL (CSV) ==================
Export.table.toDrive({
  collection: finalStats,
  description: 'Angren_NDVI_Landsat_S2_ERA5_2019_2024_seasons',
  fileFormat: 'CSV'
});

// ================== 5. КОРРЕЛЯЦИИ ==================

// 5.1. Корреляция NDVI Sentinel-2 и влажности почвы (swvl1) по всем годам и сезонам
var fc_s2_swvl1 = finalStats.filter(
  ee.Filter.notNull(['ndvi_mean_s2', 'swvl1'])
);
var corr_s2_swvl1 = fc_s2_swvl1.reduceColumns({
  reducer: ee.Reducer.pearsonsCorrelation(),
  selectors: ['ndvi_mean_s2', 'swvl1']
});
print('r (NDVI_mean_S2 vs swvl1), 2019–2024 все сезоны:', corr_s2_swvl1);

// 5.2. Корреляция NDVI Landsat и осадков (tp_mm)
var fc_ls_tp = finalStats.filter(
  ee.Filter.notNull(['ndvi_mean_ls', 'tp_mm'])
);
var corr_ls_tp = fc_ls_tp.reduceColumns({
  reducer: ee.Reducer.pearsonsCorrelation(),
  selectors: ['ndvi_mean_ls', 'tp_mm']
});
print('r (NDVI_mean_Landsat vs tp_mm), 2019–2024 все сезоны:', corr_ls_tp);

// 5.3. Пример: корреляция только для 2024 года (3 сезона)
var final2024 = finalStats.filter(ee.Filter.eq('year', 2024));

var corr2024_s2_swvl1 = final2024.filter(
  ee.Filter.notNull(['ndvi_mean_s2', 'swvl1'])
).reduceColumns({
  reducer: ee.Reducer.pearsonsCorrelation(),
  selectors: ['ndvi_mean_s2', 'swvl1']
});
print('r (NDVI_mean_S2 vs swvl1), только 2024:', corr2024_s2_swvl1);

var corr2024_ls_tp = final2024.filter(
  ee.Filter.notNull(['ndvi_mean_ls', 'tp_mm'])
).reduceColumns({
  reducer: ee.Reducer.pearsonsCorrelation(),
  selectors: ['ndvi_mean_ls', 'tp_mm']
});
print('r (NDVI_mean_Landsat vs tp_mm), только 2024:', corr2024_ls_tp); 

// ================== 6. ГРАФИКИ КОРРЕЛЯЦИЙ И ВРЕМЕННЫХ РЯДОВ ==================

// 6.1. Scatter NDVI_mean_S2 vs swvl1 (все годы и сезоны, раскраска по сезонам)
var chart_s2_swvl1_all = ui.Chart.feature.groups(
  fc_s2_swvl1,
  'ndvi_mean_s2',   // x
  'swvl1',          // y
  'season'          // группировка по сезону
).setChartType('ScatterChart')
 .setOptions({
   title: 'NDVI_mean Sentinel-2 vs swvl1 (ERA5-Land), 2019–2024, все сезоны',
   hAxis: {title: 'NDVI_mean Sentinel-2'},
   vAxis: {title: 'swvl1 (м³/м³)'},
   pointSize: 6,
   legend: {position: 'right'}
 });
print(chart_s2_swvl1_all);

// 6.2. Scatter NDVI_mean_Landsat vs tp_mm (все годы и сезоны)
var chart_ls_tp_all = ui.Chart.feature.groups(
  fc_ls_tp,
  'ndvi_mean_ls',
  'tp_mm',
  'season'
).setChartType('ScatterChart')
 .setOptions({
   title: 'NDVI_mean Landsat vs сезонные осадки (ERA5-Land), 2019–2024',
   hAxis: {title: 'NDVI_mean Landsat'},
   vAxis: {title: 'Сезонная сумма осадков, мм'},
   pointSize: 6,
   legend: {position: 'right'}
 });
print(chart_ls_tp_all);

// 6.3. Scatter только для 2024 г. (3 сезона)

// NDVI S2 vs swvl1
var final2024_s2_swvl1 = final2024.filter(
  ee.Filter.notNull(['ndvi_mean_s2', 'swvl1'])
);
var chart_s2_swvl1_2024 = ui.Chart.feature.byFeature(
  final2024_s2_swvl1,
  'ndvi_mean_s2',
  ['swvl1']
).setChartType('ScatterChart')
 .setOptions({
   title: 'NDVI_mean Sentinel-2 vs swvl1, 2024 г.',
   hAxis: {title: 'NDVI_mean Sentinel-2'},
   vAxis: {title: 'swvl1 (м³/м³)'},
   pointSize: 8
 });
print(chart_s2_swvl1_2024);

// NDVI Landsat vs tp_mm
var final2024_ls_tp = final2024.filter(
  ee.Filter.notNull(['ndvi_mean_ls', 'tp_mm'])
);
var chart_ls_tp_2024 = ui.Chart.feature.byFeature(
  final2024_ls_tp,
  'ndvi_mean_ls',
  ['tp_mm']
).setChartType('ScatterChart')
 .setOptions({
   title: 'NDVI_mean Landsat vs tp_mm, 2024 г.',
   hAxis: {title: 'NDVI_mean Landsat'},
   vAxis: {title: 'Сезонная сумма осадков, мм'},
   pointSize: 8
 });
print(chart_ls_tp_2024);

// 6.4. Временные ряды по сезонам: NDVI (Landsat + S2) и климат (tp_mm + swvl1)

var seasonNames = ['spring', 'summer', 'autumn'];

seasonNames.forEach(function(season) {
  // Фильтр по сезону
  var fcSeason = finalStats.filter(ee.Filter.eq('season', season));

  // NDVI по годам: Landsat и Sentinel-2
  var chart_ndvi_season = ui.Chart.feature.byFeature(
    fcSeason,
    'year',
    ['ndvi_mean_ls', 'ndvi_mean_s2']
  ).setChartType('LineChart')
   .setOptions({
     title: 'NDVI_mean по годам (' + season + '), Landsat и Sentinel-2',
     hAxis: {title: 'Год'},
     vAxis: {title: 'NDVI_mean'},
     lineWidth: 2,
     pointSize: 4,
     legend: {position: 'bottom'}
   });
  print(chart_ndvi_season);

  // Климат по годам: осадки и влажность почвы
  var chart_climate_season = ui.Chart.feature.byFeature(
    fcSeason,
    'year',
    ['tp_mm', 'swvl1']
  ).setChartType('LineChart')
   .setOptions({
     title: 'ERA5-Land по годам (' + season + '): осадки и влажность почвы',
     hAxis: {title: 'Год'},
     vAxis: {title: 'tp_mm (мм) / swvl1 (м³/м³)'},
     lineWidth: 2,
     pointSize: 4,
     legend: {position: 'bottom'}
   });
  print(chart_climate_season);
});
// ================== 7. ЭКСПОРТ ДАННЫХ ДЛЯ КОРРЕЛЯЦИЙ И ГРАФИКОВ ==================

// 7.1. Все годы и сезоны: NDVI_mean_S2 vs swvl1
Export.table.toDrive({
  collection: fc_s2_swvl1,
  description: 'Angren_corr_NDVImean_S2_vs_swvl1_2019_2024',
  fileFormat: 'CSV'
});

// 7.2. Все годы и сезоны: NDVI_mean_Landsat vs tp_mm
Export.table.toDrive({
  collection: fc_ls_tp,
  description: 'Angren_corr_NDVImean_Landsat_vs_tpmm_2019_2024',
  fileFormat: 'CSV'
});

// 7.3. Только 2024 год (3 сезона): полный набор NDVI + ERA5
Export.table.toDrive({
  collection: final2024,
  description: 'Angren_NDVI_Landsat_S2_ERA5_2024_seasons',
  fileFormat: 'CSV'
});

// 7.4. 2024 год: пары для scatter-графиков
Export.table.toDrive({
  collection: final2024_s2_swvl1,
  description: 'Angren_corr_2024_NDVImean_S2_vs_swvl1',
  fileFormat: 'CSV'
});

Export.table.toDrive({
  collection: final2024_ls_tp,
  description: 'Angren_corr_2024_NDVImean_Landsat_vs_tpmm',
  fileFormat: 'CSV'
});

// ================== 7. МАССОВЫЙ ЭКСПОРТ СЕЗОННЫХ КОМПОЗИТОВ NDVI ==================

// Вспомогательная функция: сезонный композит NDVI для заданного сенсора
function getSeasonNdviImage(collection, year, seasonObj) {
  year = ee.Number(year);
  seasonObj = ee.Dictionary(seasonObj);

  var startM = ee.Number(seasonObj.get('start'));
  var endM   = ee.Number(seasonObj.get('end'));

  var start = ee.Date.fromYMD(year, startM, 1);
  var end   = ee.Date.fromYMD(year, endM, 1).advance(1, 'month');

  var collSeason = collection
    .filterDate(start, end)
    .filterBounds(roi);

  // Если сцен нет – вернём null (чтобы не падал reduce/экспорт)
  var count = collSeason.size();
  print('Количество сцен', seasonObj.get('name'), year, '=>', count);

  // Средний NDVI за сезон
  var ndviSeason = collSeason.select('NDVI').mean().rename('NDVI');

  // Обрезаем к прямоугольнику 6×6 км
  return ndviSeason.clip(roi);
}

// Создаём экспорт для каждого года, сезона и сенсора
years.getInfo().forEach(function(y) {
  seasons.getInfo().forEach(function(s) {
    var seasonName = s.name;      // 'spring' / 'summer' / 'autumn'

    // -------- Landsat (30 м) --------
    var startL = ee.Date.fromYMD(y, s.start, 1);
    var endL   = ee.Date.fromYMD(y, s.end, 1).advance(1, 'month');
    var lsCollSeason = landsatAll
      .filterDate(startL, endL)
      .filterBounds(roi);

    if (lsCollSeason.size().getInfo() > 0) {
      var ndviLs = lsCollSeason.select('NDVI')
        .mean()
        .rename('NDVI')
        .clip(roi);

      Export.image.toDrive({
        image: ndviLs,
        description: 'Angren_NDVI_Landsat_' + seasonName + '_' + y + '_square6km',
        // при желании можно указать папку на Google Drive:
        // folder: 'Angren_NDVI',
        region: roi,
        scale: 30,
        maxPixels: 1e13
      });
    }

    // -------- Sentinel-2 (10 м) --------
    var startS = ee.Date.fromYMD(y, s.start, 1);
    var endS   = ee.Date.fromYMD(y, s.end, 1).advance(1, 'month');
    var s2CollSeason = S2
      .filterDate(startS, endS)
      .filterBounds(roi);

    if (s2CollSeason.size().getInfo() > 0) {
      var ndviS2 = s2CollSeason.select('NDVI')
        .mean()
        .rename('NDVI')
        .clip(roi);

      Export.image.toDrive({
        image: ndviS2,
        description: 'Angren_NDVI_S2_' + seasonName + '_' + y + '_square6km',
        // folder: 'Angren_NDVI',
        region: roi,
        scale: 10,
        maxPixels: 1e13
      });
    }
  });
});
