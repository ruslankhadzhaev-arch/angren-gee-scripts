![CodeRabbit Pull Request Reviews](https://img.shields.io/coderabbit/prs/github/ruslankhadzhaev-arch/angren-gee-scripts?utm_source=oss&utm_medium=github&utm_campaign=ruslankhadzhaev-arch%2Fangren-gee-scripts&labelColor=171717&color=FF570A&link=https%3A%2F%2Fcoderabbit.ai&label=CodeRabbit+Reviews)

/*
 * ==========================================================================================
 * ANGREN: QC + PROOF EDITION (1990–2024) — STABLE GEE JS
 * ==========================================================================================
 * BASE (QC):
 *  - MIN_PX per scene
 *  - Scene coverage mean (Cov_Scene_Mean)
 *  - Composite coverage (Cov_Composite) using pixelArea in locked projection
 * PROOF (only when QC is good):
 *  - LST_Grad_100m  : linearFit(LST ~ elevation) * 100
 *  - LST_North/South + dLST_South_North (northness masks)
 *  - NDVI by elevation zones (33/66 percentiles)
 * Climate:
 *  - ERA5-Land aggregated at native scale ~11 km (no resampling to 30 m)
 * ==========================================================================================
 */

// -----------------------------
// 1) SETTINGS
// -----------------------------
var roi = ee.Geometry.Rectangle([70.1446, 41.0091, 70.2158, 41.0632]);
Map.centerObject(roi, 12);
Map.addLayer(roi, {color: 'red', fillColor: '00000000'}, 'ROI Angren');

var startYear = 1990;
var endYear   = 2024;
var years     = ee.List.sequence(startYear, endYear);

var seasons = [
  {name: 'Spring', start: 3, end: 5,  order: 1},
  {name: 'Summer', start: 6, end: 8,  order: 2},
  {name: 'Autumn', start: 9, end: 11, order: 3}
];

var MIN_PX = 10;                 // minimum valid pixels per scene
var MIN_SCENES_PROOF = 3;        // minimum scenes to compute proof metrics
var COV_THRESHOLD = 0.5;         // minimum composite coverage to compute proof metrics

// -----------------------------
// 2) REFERENCE PROJECTION (LOCK) + ROI AREA IN THAT PROJECTION
// -----------------------------
var refImg = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
  .merge(ee.ImageCollection('LANDSAT/LE07/C02/T1_L2'))
  .merge(ee.ImageCollection('LANDSAT/LC08/C02/T1_L2'))
  .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2'))
  .filterBounds(roi)
  .first();

// If no scenes exist in ROI (rare), this might fail without checks, but assumed stable for this ROI.
var refProj = ee.Image(refImg).select('QA_PIXEL').projection();

var roiAreaGeo = ee.Number(roi.area(1));
var roiArea30m = ee.Number(
  ee.Image.pixelArea().reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: roi,
    scale: 30,
    crs: refProj,
    maxPixels: 1e9,
    tileScale: 4
  }).get('area')
);

// -----------------------------
// 3) HELPERS: Modern sensor check & Masking
// -----------------------------
// Robust check: Landsat 8/9 have thermal band ST_B10, L5/7 do not.
function isL8L9(image) {
  return image.bandNames().contains('ST_B10');
}

/**
 * Strict QA_PIXEL:
 * - Always exclude bits: 0 Fill, 1 Dilated Cloud, 3 Cloud, 4 Shadow, 5 Snow/Ice
 * - Exclude bit 2 Cirrus ONLY for Landsat 8/9
 */
function maskLandsatStrict(image) {
  var qa = image.select('QA_PIXEL');

  var baseMask = qa.bitwiseAnd(1 << 0).eq(0)
    .and(qa.bitwiseAnd(1 << 1).eq(0))
    .and(qa.bitwiseAnd(1 << 3).eq(0))
    .and(qa.bitwiseAnd(1 << 4).eq(0))
    .and(qa.bitwiseAnd(1 << 5).eq(0));

  var modern = isL8L9(image);
  var finalMask = ee.Image(ee.Algorithms.If(
    modern,
    baseMask.and(qa.bitwiseAnd(1 << 2).eq(0)), // cirrus only for L8/9
    baseMask
  ));

  return image.updateMask(finalMask);
}

// -----------------------------
// 4) SCALE FACTORS + BAND HARMONIZATION + INDICES
// -----------------------------
function applyScaleFactors(image) {
  var optical = image.select('SR_B.').multiply(0.0000275).add(-0.2);
  var thermal = image.select('ST_B.*').multiply(0.00341802).add(149.0);
  return image.addBands(optical, null, true)
              .addBands(thermal, null, true);
}

function unifyBands(image) {
  var modern = isL8L9(image);

  var bands = ee.List(ee.Algorithms.If(
    modern,
    ['SR_B4', 'SR_B5', 'ST_B10'],  // L8/9
    ['SR_B3', 'SR_B4', 'ST_B6']    // L5/7
  ));

  return image.select(bands, ['Red', 'NIR', 'Temp'])
    .copyProperties(image, ['system:time_start', 'LANDSAT_SCENE_ID', 'SPACECRAFT_ID']);
}

function addIndices(image) {
  var ndvi = image.normalizedDifference(['NIR', 'Red']).rename('NDVI');
  var savi = image.expression(
    '((NIR - RED) / (NIR + RED + 0.5)) * 1.5',
    {'NIR': image.select('NIR'), 'RED': image.select('Red')}
  ).rename('SAVI');

  // Temp (scaled) is Kelvin -> Celsius
  var lstC = image.select('Temp').subtract(273.15).rename('LST_C');

  return image.addBands([ndvi, savi, lstC]);
}

// -----------------------------
// 5) TERRAIN (NASADEM) + ZONES
// -----------------------------
var dem = ee.Image('NASA/NASADEM_HGT/001');
var elevation = dem.select('elevation');
var slope = ee.Terrain.slope(elevation).rename('Slope');

var aspectRad = ee.Terrain.aspect(elevation).multiply(Math.PI).divide(180);
var northness = aspectRad.cos().rename('Northness'); // +1 north, -1 south

var terrainStack = elevation.rename('elevation')
  .addBands(slope)
  .addBands(northness);

// Elevation percentiles for zones (computed once for ROI)
var elevPerc = elevation.reduceRegion({
  reducer: ee.Reducer.percentile([33, 66]),
  geometry: roi,
  scale: 30,
  crs: refProj,
  maxPixels: 1e9,
  tileScale: 4
});
var p33 = ee.Number(elevPerc.get('elevation_p33'));
var p66 = ee.Number(elevPerc.get('elevation_p66'));

var zoneLow  = elevation.lt(p33);
var zoneMid  = elevation.gte(p33).and(elevation.lt(p66));
var zoneHigh = elevation.gte(p66);

// -----------------------------
// 6) ERA5-LAND (native ~11 km) aggregation per year-season
// -----------------------------
function aggregateClimate(year, season) {
  var eraCol = ee.ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR')
    .filterBounds(roi)
    .filter(ee.Filter.calendarRange(year, year, 'year'))
    .filter(ee.Filter.calendarRange(season.start, season.end, 'month'));

  // precipitation: seasonal sum (m) -> spatial mean over ROI
  var pp_m = eraCol.select('total_precipitation_sum').sum()
    .reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: roi,
      scale: 11132,
      maxPixels: 1e9
    }).get('total_precipitation_sum');

  pp_m = ee.Algorithms.If(pp_m, pp_m, 0);
  var pp_mm = ee.Number(pp_m).multiply(1000);

  // soil moisture: seasonal mean -> spatial mean over ROI
  var ssm = eraCol.select('volumetric_soil_water_layer_1').mean()
    .reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: roi,
      scale: 11132,
      maxPixels: 1e9
    }).get('volumetric_soil_water_layer_1');

  ssm = ee.Algorithms.If(ssm, ssm, 0);

  return ee.Dictionary({'PP_mm': pp_mm, 'SSM_m3m3': ssm});
}

// -----------------------------
// 7) HELPERS for robust reduceRegion means with null
// -----------------------------
function meanBand(img, bandName, maskImg) {
  var x = img.select(bandName);
  x = ee.Image(ee.Algorithms.If(maskImg, x.updateMask(maskImg), x));

  var v = x.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: roi,
    scale: 30,
    crs: refProj,
    maxPixels: 1e9,
    tileScale: 4
  }).get(bandName);

  return ee.Algorithms.If(v, v, null);
}

// -----------------------------
// 8) MAIN LOOP: YEAR x SEASON
// -----------------------------
var finalResults = ee.FeatureCollection(
  years.map(function(y) {
    return ee.List(seasons).map(function(seasonObj) {
      var season = ee.Dictionary(seasonObj);
      var year = ee.Number(y);

      // Climate
      var climate = aggregateClimate(year, {
        start: season.getNumber('start'),
        end: season.getNumber('end')
      });

      // Landsat collection
      var col = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
        .merge(ee.ImageCollection('LANDSAT/LE07/C02/T1_L2'))
        .merge(ee.ImageCollection('LANDSAT/LC08/C02/T1_L2'))
        .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2'))
        .filterBounds(roi)
        .filter(ee.Filter.calendarRange(year, year, 'year'))
        .filter(ee.Filter.calendarRange(season.getNumber('start'), season.getNumber('end'), 'month'));

      // Step 1: valid pixel count after strict mask (on QA mask)
      var withPx = col.map(function(img) {
        var masked = maskLandsatStrict(img);
        var px = masked.select('QA_PIXEL').reduceRegion({
          reducer: ee.Reducer.count(),
          geometry: roi,
          scale: 30,
          crs: refProj,
          maxPixels: 1e9,
          tileScale: 4
        }).get('QA_PIXEL');
        px = ee.Number(ee.Algorithms.If(px, px, 0));
        return img.set('valid_pixels', px);
      });

      var validImages = withPx.filter(ee.Filter.gte('valid_pixels', MIN_PX));
      var count = validImages.size();

      // Step 2: per-scene coverage ratio
      var validWithCoverage = validImages.map(function(img) {
        var masked = maskLandsatStrict(img);
        var sceneArea = ee.Image.pixelArea()
          .updateMask(masked.select('QA_PIXEL').mask())
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

      var covSceneMean = ee.Number(ee.Algorithms.If(
        count.gt(0),
        validWithCoverage.aggregate_mean('scene_cov_area_ratio'),
        0
      ));

      // If count == 0, return an empty row with climate
      return ee.Algorithms.If(count.gt(0), (function() {

        // processed collection -> composite
        var processed = ee.ImageCollection(validWithCoverage)
          .map(maskLandsatStrict)
          .map(applyScaleFactors)
          .map(unifyBands)
          .map(addIndices);

        var composite = processed.median().clip(roi);

        // composite coverage (use Red mask)
        var compArea = ee.Image.pixelArea()
          .updateMask(composite.select('Red').mask())
          .reduceRegion({
            reducer: ee.Reducer.sum(),
            geometry: roi,
            scale: 30,
            crs: refProj,
            maxPixels: 1e9,
            tileScale: 4
          }).get('area');

        compArea = ee.Number(ee.Algorithms.If(compArea, compArea, 0));
        var covComp = compArea.divide(roiArea30m).min(1).max(0);

        // Stats for indices
        var stats = composite.select(['NDVI','SAVI','LST_C']).reduceRegion({
          reducer: ee.Reducer.mean().combine({
            reducer2: ee.Reducer.stdDev(),
            sharedInputs: true
          }),
          geometry: roi,
          scale: 30,
          crs: refProj,
          maxPixels: 1e9,
          tileScale: 4
        });

        // QC condition for PROOF metrics
        var cond = covComp.gte(COV_THRESHOLD).and(count.gte(MIN_SCENES_PROOF));

        // analysis stack for proof
        var analysisStack = composite.addBands(terrainStack);

        // PROOF 1: vertical thermal gradient (°C/100m)
        // Fixed: Wrapped in conditional to match "PROOF (only when QC is good)"
        var lstGrad100m = ee.Algorithms.If(cond, (function() {
           var regression = analysisStack.select(['elevation', 'LST_C']).reduceRegion({
             reducer: ee.Reducer.linearFit(),
             geometry: roi,
             scale: 30,
             crs: refProj,
             maxPixels: 1e9,
             tileScale: 4
           });
           return ee.Number(ee.Algorithms.If(
             regression.get('scale'),
             ee.Number(regression.get('scale')).multiply(100),
             null
           ));
        })(), null);

        // PROOF 2: South vs North slopes (optionally require slope > 5°)
        var maskNorth = northness.gt(0.3); // .and(slope.gt(5));
        var maskSouth = northness.lt(-0.3); // .and(slope.gt(5));

        var LST_North = ee.Algorithms.If(
          cond,
          meanBand(analysisStack, 'LST_C', maskNorth),
          null
        );

        var LST_South = ee.Algorithms.If(
          cond,
          meanBand(analysisStack, 'LST_C', maskSouth),
          null
        );

        // dLST considered only if cond=true AND both values really exist
        var dLST_South_North = ee.Algorithms.If(
          cond,
          ee.Algorithms.If(
            LST_North,
            ee.Algorithms.If(
              LST_South,
              ee.Number(LST_South).subtract(ee.Number(LST_North)),
              null
            ),
            null
          ),
          null
        );

        // PROOF 3: NDVI by elevation zones
        var NDVI_Low_Zone = ee.Algorithms.If(cond, meanBand(analysisStack, 'NDVI', zoneLow), null);
        var NDVI_Mid_Zone = ee.Algorithms.If(cond, meanBand(analysisStack, 'NDVI', zoneMid), null);
        var NDVI_High_Zone = ee.Algorithms.If(cond, meanBand(analysisStack, 'NDVI', zoneHigh), null);

        return ee.Feature(null, {
          'Year': year,
          'Season': season.getString('name'),
          'SeasonOrder': season.getNumber('order'),

          'Image_Count': count,
          'MIN_PX': MIN_PX,
          'ROI_Area_Geo_m2': roiAreaGeo,
          'ROI_Area_30m_m2': roiArea30m,

          'Cov_Scene_Mean': covSceneMean,
          'Cov_Composite': covComp,

          'NDVI_mean': stats.get('NDVI_mean'),
          'NDVI_std': stats.get('NDVI_stdDev'),
          'SAVI_mean': stats.get('SAVI_mean'),
          'SAVI_std': stats.get('SAVI_stdDev'),
          'LST_mean_C': stats.get('LST_C_mean'),
          'LST_std_C': stats.get('LST_C_stdDev'),

          'PP_mm': climate.get('PP_mm'),
          'SSM_m3m3': climate.get('SSM_m3m3'),

          // PROOF metrics
          'LST_Grad_100m': lstGrad100m,
          'LST_North': LST_North,
          'LST_South': LST_South,
          'dLST_South_North': dLST_South_North,
          'NDVI_Low_Zone': NDVI_Low_Zone,
          'NDVI_Mid_Zone': NDVI_Mid_Zone,
          'NDVI_High_Zone': NDVI_High_Zone
        });

      })(), ee.Feature(null, {
        'Year': year,
        'Season': season.getString('name'),
        'SeasonOrder': season.getNumber('order'),
        'Image_Count': 0,
        'MIN_PX': MIN_PX,
        'Cov_Scene_Mean': 0,
        'Cov_Composite': 0,
        'PP_mm': climate.get('PP_mm'),
        'SSM_m3m3': climate.get('SSM_m3m3')
      }));

    });
  }).flatten()
);

// Sort
var finalData = finalResults.sort('Year').sort('SeasonOrder');

// -----------------------------
// 9) CHARTS (examples)
// -----------------------------
print('QC + PROOF dataset:', finalData.limit(10));

var summer = finalData.filter(ee.Filter.eq('Season', 'Summer'));

print(ui.Chart.feature.byFeature(summer, 'Year', ['NDVI_mean'])
  .setOptions({
    title: 'NDVI_mean (Summer)',
    vAxis: {title: 'NDVI'},
    interpolateNulls: false,
    pointSize: 4
  })
);

print(ui.Chart.feature.byFeature(summer, 'Year', ['LST_Grad_100m'])
  .setOptions({
    title: 'LST vertical gradient (Summer) — °C per 100 m (QC-gated)',
    vAxis: {title: '°C per 100 m'},
    interpolateNulls: false,
    pointSize: 4
  })
);

// -----------------------------
// 10) EXPORT
// -----------------------------
Export.table.toDrive({
  collection: finalData,
  description: 'Angren_QC_Proof_Dataset',
  fileFormat: 'CSV',
  selectors: [
    'Year','Season','SeasonOrder',
    'Image_Count','MIN_PX',
    'ROI_Area_Geo_m2','ROI_Area_30m_m2',
    'Cov_Scene_Mean','Cov_Composite',
    'NDVI_mean','NDVI_std',
    'SAVI_mean','SAVI_std',
    'LST_mean_C','LST_std_C',
    'PP_mm','SSM_m3m3',
    'LST_Grad_100m','LST_North','LST_South','dLST_South_North',
    'NDVI_Low_Zone','NDVI_Mid_Zone','NDVI_High_Zone'
  ]
});

// =====================================================================
// FULL EXPORT (PER YEAR × PER SEASON) + STUDY AREA PACK — FIXED CRS
// Key fix: DO NOT pass ee.Projection (refProj) into Export "crs".
// Use a CRS STRING (recommended) or omit "crs" entirely.
// =====================================================================

// ----------------------------
// 0) EXPORT SETTINGS (edit here)
// ----------------------------
var EXPORT_FOLDER = 'Angren_JournalPack';
var MAP_SCALE = 30;

// ✅ Recommended fixed CRS for Angren (~70E, 41N) = UTM 42N
// If you prefer to keep native image projection, set EXPORT_CRS = null
var EXPORT_CRS = 'EPSG:32642';   // or set to null to omit crs

var MAPS_ONLY_WHEN_DATA_EXISTS = true;

var MAPS_QC_GATE = true;
var MAPS_MIN_COV = 0.80;
var MAPS_MIN_SCENES = 1;

var EXPORT_NDVI = true;
var EXPORT_LST  = true;
var EXPORT_MASK = true;

var EXPORT_STUDYAREA_VECTOR = true;
var EXPORT_STUDYAREA_HILLSHADE = true;
var EXPORT_STUDYAREA_DEM = true;
var EXPORT_STUDYAREA_STATIC_TERRAIN_ZONES = true;

var seasonMonths = {
  'Spring': {start: 3, end: 5},
  'Summer': {start: 6, end: 8},
  'Autumn': {start: 9, end: 11}
};

// ----------------------------
// 1) STUDY AREA PACK (Figure 1 base)
// ----------------------------
if (EXPORT_STUDYAREA_VECTOR) {
  var roiFc = ee.FeatureCollection([ee.Feature(roi, {'name': 'ROI_Angren'})]);
  Export.table.toDrive({
    collection: roiFc,
    description: 'Angren_StudyArea_ROI',
    folder: EXPORT_FOLDER,
    fileFormat: 'SHP'
  });
}

if (EXPORT_STUDYAREA_HILLSHADE) {
  var hillshade = ee.Terrain.hillshade(elevation).rename('Hillshade');

  var hsParams = {
    image: hillshade.toUint8(),
    description: 'Angren_StudyArea_Hillshade',
    folder: EXPORT_FOLDER,
    scale: MAP_SCALE,
    region: roi,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF'
  };
  if (EXPORT_CRS) hsParams.crs = EXPORT_CRS;

  Export.image.toDrive(hsParams);
}

if (EXPORT_STUDYAREA_DEM) {
  var demParams = {
    image: elevation.rename('Elevation_m').toFloat(),
    description: 'Angren_StudyArea_DEM',
    folder: EXPORT_FOLDER,
    scale: MAP_SCALE,
    region: roi,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF'
  };
  if (EXPORT_CRS) demParams.crs = EXPORT_CRS;

  Export.image.toDrive(demParams);
}

if (EXPORT_STUDYAREA_STATIC_TERRAIN_ZONES) {
  var staticMaps = terrainStack
    .addBands(zoneLow.rename('Zone_Low').toUint8())
    .addBands(zoneMid.rename('Zone_Mid').toUint8())
    .addBands(zoneHigh.rename('Zone_High').toUint8());

  var stParams = {
    image: staticMaps,
    description: 'Angren_Static_Terrain_Zones',
    folder: EXPORT_FOLDER,
    scale: MAP_SCALE,
    region: roi,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF'
  };
  if (EXPORT_CRS) stParams.crs = EXPORT_CRS;

  Export.image.toDrive(stParams);
}

// ----------------------------
// 2) FULL MAP EXPORTS: each Year × Season (NDVI, LST_C, ValidMask)
// ----------------------------
function makeComposite(year, seasonName) {
  var m = seasonMonths[seasonName];

  var col = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
    .merge(ee.ImageCollection('LANDSAT/LE07/C02/T1_L2'))
    .merge(ee.ImageCollection('LANDSAT/LC08/C02/T1_L2'))
    .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2'))
    .filterBounds(roi)
    .filter(ee.Filter.calendarRange(year, year, 'year'))
    .filter(ee.Filter.calendarRange(m.start, m.end, 'month'))
    .map(maskLandsatStrict)
    .map(applyScaleFactors)
    .map(unifyBands)
    .map(addIndices);

  return col.median().clip(roi);
}

function exportAllYearSeasonMapsFromFinalData() {
  var fc = finalData;

  if (MAPS_ONLY_WHEN_DATA_EXISTS) {
    fc = fc.filter(ee.Filter.gte('Image_Count', MAPS_MIN_SCENES));
  }
  if (MAPS_QC_GATE) {
    fc = fc.filter(ee.Filter.gte('Cov_Composite', MAPS_MIN_COV));
  }

  fc = fc.select(['Year', 'Season', 'Image_Count', 'Cov_Composite']);

  fc.toList(fc.size()).evaluate(function(list) {
    if (!list || list.length === 0) {
      print('No year-season combinations passed the export filters.');
      return;
    }

    print('Year-season exports to be created:', list.length);

    list.forEach(function(f) {
      var p = f.properties;
      var year = Number(p.Year);
      var seasonName = String(p.Season);

      if (!seasonMonths[seasonName]) return;

      var comp = makeComposite(year, seasonName);
      var tag = seasonName + '_' + year;

      if (EXPORT_NDVI) {
        var ndviParams = {
          image: comp.select('NDVI').toFloat(),
          description: 'Angren_NDVI_' + tag,
          folder: EXPORT_FOLDER,
          scale: MAP_SCALE,
          region: roi,
          maxPixels: 1e13,
          fileFormat: 'GeoTIFF'
        };
        if (EXPORT_CRS) ndviParams.crs = EXPORT_CRS;
        Export.image.toDrive(ndviParams);
      }

      if (EXPORT_LST) {
        var lstParams = {
          image: comp.select('LST_C').toFloat(),
          description: 'Angren_LST_' + tag,
          folder: EXPORT_FOLDER,
          scale: MAP_SCALE,
          region: roi,
          maxPixels: 1e13,
          fileFormat: 'GeoTIFF'
        };
        if (EXPORT_CRS) lstParams.crs = EXPORT_CRS;
        Export.image.toDrive(lstParams);
      }

      if (EXPORT_MASK) {
        var maskParams = {
          image: comp.select('Red').mask().rename('ValidMask').toUint8(),
          description: 'Angren_ValidMask_' + tag,
          folder: EXPORT_FOLDER,
          scale: MAP_SCALE,
          region: roi,
          maxPixels: 1e13,
          fileFormat: 'GeoTIFF'
        };
        if (EXPORT_CRS) maskParams.crs = EXPORT_CRS;
        Export.image.toDrive(maskParams);
      }
    });

    print('All export tasks created. Open the Tasks tab and run them.');
  });
}

// =====================================================================
// FULL EXPORT (PER YEAR × PER SEASON) + STUDY AREA PACK — FIXED CRS
// Key fix: DO NOT pass ee.Projection (refProj) into Export "crs".
// Use a CRS STRING (recommended) or omit "crs" entirely.
// =====================================================================

// ----------------------------
// 0) EXPORT SETTINGS (edit here)
// ----------------------------
var EXPORT_FOLDER = 'Angren_JournalPack';
var MAP_SCALE = 30;

// ✅ Recommended fixed CRS for Angren (~70E, 41N) = UTM 42N
// If you prefer to keep native image projection, set EXPORT_CRS = null
var EXPORT_CRS = 'EPSG:32642';   // or set to null to omit crs

var MAPS_ONLY_WHEN_DATA_EXISTS = true;

var MAPS_QC_GATE = true;
var MAPS_MIN_COV = 0.80;
var MAPS_MIN_SCENES = 1;

var EXPORT_NDVI = true;
var EXPORT_LST  = true;
var EXPORT_MASK = true;

var EXPORT_STUDYAREA_VECTOR = true;
var EXPORT_STUDYAREA_HILLSHADE = true;
var EXPORT_STUDYAREA_DEM = true;
var EXPORT_STUDYAREA_STATIC_TERRAIN_ZONES = true;

var seasonMonths = {
  'Spring': {start: 3, end: 5},
  'Summer': {start: 6, end: 8},
  'Autumn': {start: 9, end: 11}
};

// ----------------------------
// 1) STUDY AREA PACK (Figure 1 base)
// ----------------------------
if (EXPORT_STUDYAREA_VECTOR) {
  var roiFc = ee.FeatureCollection([ee.Feature(roi, {'name': 'ROI_Angren'})]);
  Export.table.toDrive({
    collection: roiFc,
    description: 'Angren_StudyArea_ROI',
    folder: EXPORT_FOLDER,
    fileFormat: 'SHP'
  });
}

if (EXPORT_STUDYAREA_HILLSHADE) {
  var hillshade = ee.Terrain.hillshade(elevation).rename('Hillshade');

  var hsParams = {
    image: hillshade.toUint8(),
    description: 'Angren_StudyArea_Hillshade',
    folder: EXPORT_FOLDER,
    scale: MAP_SCALE,
    region: roi,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF'
  };
  if (EXPORT_CRS) hsParams.crs = EXPORT_CRS;

  Export.image.toDrive(hsParams);
}

if (EXPORT_STUDYAREA_DEM) {
  var demParams = {
    image: elevation.rename('Elevation_m').toFloat(),
    description: 'Angren_StudyArea_DEM',
    folder: EXPORT_FOLDER,
    scale: MAP_SCALE,
    region: roi,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF'
  };
  if (EXPORT_CRS) demParams.crs = EXPORT_CRS;

  Export.image.toDrive(demParams);
}

if (EXPORT_STUDYAREA_STATIC_TERRAIN_ZONES) {
  var staticMaps = terrainStack
    .addBands(zoneLow.rename('Zone_Low').toUint8())
    .addBands(zoneMid.rename('Zone_Mid').toUint8())
    .addBands(zoneHigh.rename('Zone_High').toUint8());

  var stParams = {
    image: staticMaps,
    description: 'Angren_Static_Terrain_Zones',
    folder: EXPORT_FOLDER,
    scale: MAP_SCALE,
    region: roi,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF'
  };
  if (EXPORT_CRS) stParams.crs = EXPORT_CRS;

  Export.image.toDrive(stParams);
}

// ----------------------------
// 2) FULL MAP EXPORTS: each Year × Season (NDVI, LST_C, ValidMask)
// ----------------------------
function makeComposite(year, seasonName) {
  var m = seasonMonths[seasonName];

  var col = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
    .merge(ee.ImageCollection('LANDSAT/LE07/C02/T1_L2'))
    .merge(ee.ImageCollection('LANDSAT/LC08/C02/T1_L2'))
    .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2'))
    .filterBounds(roi)
    .filter(ee.Filter.calendarRange(year, year, 'year'))
    .filter(ee.Filter.calendarRange(m.start, m.end, 'month'))
    .map(maskLandsatStrict)
    .map(applyScaleFactors)
    .map(unifyBands)
    .map(addIndices);

  return col.median().clip(roi);
}

function exportAllYearSeasonMapsFromFinalData() {
  var fc = finalData;

  if (MAPS_ONLY_WHEN_DATA_EXISTS) {
    fc = fc.filter(ee.Filter.gte('Image_Count', MAPS_MIN_SCENES));
  }
  if (MAPS_QC_GATE) {
    fc = fc.filter(ee.Filter.gte('Cov_Composite', MAPS_MIN_COV));
  }

  fc = fc.select(['Year', 'Season', 'Image_Count', 'Cov_Composite']);

  fc.toList(fc.size()).evaluate(function(list) {
    if (!list || list.length === 0) {
      print('No year-season combinations passed the export filters.');
      return;
    }

    print('Year-season exports to be created:', list.length);

    list.forEach(function(f) {
      var p = f.properties;
      var year = Number(p.Year);
      var seasonName = String(p.Season);

      if (!seasonMonths[seasonName]) return;

      var comp = makeComposite(year, seasonName);
      var tag = seasonName + '_' + year;

      if (EXPORT_NDVI) {
        var ndviParams = {
          image: comp.select('NDVI').toFloat(),
          description: 'Angren_NDVI_' + tag,
          folder: EXPORT_FOLDER,
          scale: MAP_SCALE,
          region: roi,
          maxPixels: 1e13,
          fileFormat: 'GeoTIFF'
        };
        if (EXPORT_CRS) ndviParams.crs = EXPORT_CRS;
        Export.image.toDrive(ndviParams);
      }

      if (EXPORT_LST) {
        var lstParams = {
          image: comp.select('LST_C').toFloat(),
          description: 'Angren_LST_' + tag,
          folder: EXPORT_FOLDER,
          scale: MAP_SCALE,
          region: roi,
          maxPixels: 1e13,
          fileFormat: 'GeoTIFF'
        };
        if (EXPORT_CRS) lstParams.crs = EXPORT_CRS;
        Export.image.toDrive(lstParams);
      }

      if (EXPORT_MASK) {
        var maskParams = {
          image: comp.select('Red').mask().rename('ValidMask').toUint8(),
          description: 'Angren_ValidMask_' + tag,
          folder: EXPORT_FOLDER,
          scale: MAP_SCALE,
          region: roi,
          maxPixels: 1e13,
          fileFormat: 'GeoTIFF'
        };
        if (EXPORT_CRS) maskParams.crs = EXPORT_CRS;
        Export.image.toDrive(maskParams);
      }
    });

    print('All export tasks created. Open the Tasks tab and run them.');
  });
}

// =====================================================================
// FULL EXPORT (PER YEAR × PER SEASON) + STUDY AREA PACK — FIXED CRS
// Key fix: DO NOT pass ee.Projection (refProj) into Export "crs".
// Use a CRS STRING (recommended) or omit "crs" entirely.
// =====================================================================

// ----------------------------
// 0) EXPORT SETTINGS (edit here)
// ----------------------------
var EXPORT_FOLDER = 'Angren_JournalPack';
var MAP_SCALE = 30;

// ✅ Recommended fixed CRS for Angren (~70E, 41N) = UTM 42N
// If you prefer to keep native image projection, set EXPORT_CRS = null
var EXPORT_CRS = 'EPSG:32642';   // or set to null to omit crs

var MAPS_ONLY_WHEN_DATA_EXISTS = true;

var MAPS_QC_GATE = true;
var MAPS_MIN_COV = 0.80;
var MAPS_MIN_SCENES = 1;

var EXPORT_NDVI = true;
var EXPORT_LST  = true;
var EXPORT_MASK = true;

var EXPORT_STUDYAREA_VECTOR = true;
var EXPORT_STUDYAREA_HILLSHADE = true;
var EXPORT_STUDYAREA_DEM = true;
var EXPORT_STUDYAREA_STATIC_TERRAIN_ZONES = true;

var seasonMonths = {
  'Spring': {start: 3, end: 5},
  'Summer': {start: 6, end: 8},
  'Autumn': {start: 9, end: 11}
};

// ----------------------------
// 1) STUDY AREA PACK (Figure 1 base)
// ----------------------------
if (EXPORT_STUDYAREA_VECTOR) {
  var roiFc = ee.FeatureCollection([ee.Feature(roi, {'name': 'ROI_Angren'})]);
  Export.table.toDrive({
    collection: roiFc,
    description: 'Angren_StudyArea_ROI',
    folder: EXPORT_FOLDER,
    fileFormat: 'SHP'
  });
}

if (EXPORT_STUDYAREA_HILLSHADE) {
  var hillshade = ee.Terrain.hillshade(elevation).rename('Hillshade');

  var hsParams = {
    image: hillshade.toUint8(),
    description: 'Angren_StudyArea_Hillshade',
    folder: EXPORT_FOLDER,
    scale: MAP_SCALE,
    region: roi,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF'
  };
  if (EXPORT_CRS) hsParams.crs = EXPORT_CRS;

  Export.image.toDrive(hsParams);
}

if (EXPORT_STUDYAREA_DEM) {
  var demParams = {
    image: elevation.rename('Elevation_m').toFloat(),
    description: 'Angren_StudyArea_DEM',
    folder: EXPORT_FOLDER,
    scale: MAP_SCALE,
    region: roi,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF'
  };
  if (EXPORT_CRS) demParams.crs = EXPORT_CRS;

  Export.image.toDrive(demParams);
}

if (EXPORT_STUDYAREA_STATIC_TERRAIN_ZONES) {
  var staticMaps = terrainStack
    .addBands(zoneLow.rename('Zone_Low').toUint8())
    .addBands(zoneMid.rename('Zone_Mid').toUint8())
    .addBands(zoneHigh.rename('Zone_High').toUint8());

  var stParams = {
    image: staticMaps,
    description: 'Angren_Static_Terrain_Zones',
    folder: EXPORT_FOLDER,
    scale: MAP_SCALE,
    region: roi,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF'
  };
  if (EXPORT_CRS) stParams.crs = EXPORT_CRS;

  Export.image.toDrive(stParams);
}

// ----------------------------
// 2) FULL MAP EXPORTS: each Year × Season (NDVI, LST_C, ValidMask)
// ----------------------------
function makeComposite(year, seasonName) {
  var m = seasonMonths[seasonName];

  var col = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
    .merge(ee.ImageCollection('LANDSAT/LE07/C02/T1_L2'))
    .merge(ee.ImageCollection('LANDSAT/LC08/C02/T1_L2'))
    .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2'))
    .filterBounds(roi)
    .filter(ee.Filter.calendarRange(year, year, 'year'))
    .filter(ee.Filter.calendarRange(m.start, m.end, 'month'))
    .map(maskLandsatStrict)
    .map(applyScaleFactors)
    .map(unifyBands)
    .map(addIndices);

  return col.median().clip(roi);
}

function exportAllYearSeasonMapsFromFinalData() {
  var fc = finalData;

  if (MAPS_ONLY_WHEN_DATA_EXISTS) {
    fc = fc.filter(ee.Filter.gte('Image_Count', MAPS_MIN_SCENES));
  }
  if (MAPS_QC_GATE) {
    fc = fc.filter(ee.Filter.gte('Cov_Composite', MAPS_MIN_COV));
  }

  fc = fc.select(['Year', 'Season', 'Image_Count', 'Cov_Composite']);

  fc.toList(fc.size()).evaluate(function(list) {
    if (!list || list.length === 0) {
      print('No year-season combinations passed the export filters.');
      return;
    }

    print('Year-season exports to be created:', list.length);

    list.forEach(function(f) {
      var p = f.properties;
      var year = Number(p.Year);
      var seasonName = String(p.Season);

      if (!seasonMonths[seasonName]) return;

      var comp = makeComposite(year, seasonName);
      var tag = seasonName + '_' + year;

      if (EXPORT_NDVI) {
        var ndviParams = {
          image: comp.select('NDVI').toFloat(),
          description: 'Angren_NDVI_' + tag,
          folder: EXPORT_FOLDER,
          scale: MAP_SCALE,
          region: roi,
          maxPixels: 1e13,
          fileFormat: 'GeoTIFF'
        };
        if (EXPORT_CRS) ndviParams.crs = EXPORT_CRS;
        Export.image.toDrive(ndviParams);
      }

      if (EXPORT_LST) {
        var lstParams = {
          image: comp.select('LST_C').toFloat(),
          description: 'Angren_LST_' + tag,
          folder: EXPORT_FOLDER,
          scale: MAP_SCALE,
          region: roi,
          maxPixels: 1e13,
          fileFormat: 'GeoTIFF'
        };
        if (EXPORT_CRS) lstParams.crs = EXPORT_CRS;
        Export.image.toDrive(lstParams);
      }

      if (EXPORT_MASK) {
        var maskParams = {
          image: comp.select('Red').mask().rename('ValidMask').toUint8(),
          description: 'Angren_ValidMask_' + tag,
          folder: EXPORT_FOLDER,
          scale: MAP_SCALE,
          region: roi,
          maxPixels: 1e13,
          fileFormat: 'GeoTIFF'
        };
        if (EXPORT_CRS) maskParams.crs = EXPORT_CRS;
        Export.image.toDrive(maskParams);
      }
    });

    print('All export tasks created. Open the Tasks tab and run them.');
  });
}

// =====================================================================
// FULL EXPORT (PER YEAR × PER SEASON) + STUDY AREA PACK — FIXED CRS
// Key fix: DO NOT pass ee.Projection (refProj) into Export "crs".
// Use a CRS STRING (recommended) or omit "crs" entirely.
// =====================================================================

// ----------------------------
// 0) EXPORT SETTINGS (edit here)
// ----------------------------
var EXPORT_FOLDER = 'Angren_JournalPack';
var MAP_SCALE = 30;

// ✅ Recommended fixed CRS for Angren (~70E, 41N) = UTM 42N
// If you prefer to keep native image projection, set EXPORT_CRS = null
var EXPORT_CRS = 'EPSG:32642';   // or set to null to omit crs

var MAPS_ONLY_WHEN_DATA_EXISTS = true;

var MAPS_QC_GATE = true;
var MAPS_MIN_COV = 0.80;
var MAPS_MIN_SCENES = 1;

var EXPORT_NDVI = true;
var EXPORT_LST  = true;
var EXPORT_MASK = true;

var EXPORT_STUDYAREA_VECTOR = true;
var EXPORT_STUDYAREA_HILLSHADE = true;
var EXPORT_STUDYAREA_DEM = true;
var EXPORT_STUDYAREA_STATIC_TERRAIN_ZONES = true;

var seasonMonths = {
  'Spring': {start: 3, end: 5},
  'Summer': {start: 6, end: 8},
  'Autumn': {start: 9, end: 11}
};

// ----------------------------
// 1) STUDY AREA PACK (Figure 1 base)
// ----------------------------
if (EXPORT_STUDYAREA_VECTOR) {
  var roiFc = ee.FeatureCollection([ee.Feature(roi, {'name': 'ROI_Angren'})]);
  Export.table.toDrive({
    collection: roiFc,
    description: 'Angren_StudyArea_ROI',
    folder: EXPORT_FOLDER,
    fileFormat: 'SHP'
  });
}

if (EXPORT_STUDYAREA_HILLSHADE) {
  var hillshade = ee.Terrain.hillshade(elevation).rename('Hillshade');

  var hsParams = {
    image: hillshade.toUint8(),
    description: 'Angren_StudyArea_Hillshade',
    folder: EXPORT_FOLDER,
    scale: MAP_SCALE,
    region: roi,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF'
  };
  if (EXPORT_CRS) hsParams.crs = EXPORT_CRS;

  Export.image.toDrive(hsParams);
}

if (EXPORT_STUDYAREA_DEM) {
  var demParams = {
    image: elevation.rename('Elevation_m').toFloat(),
    description: 'Angren_StudyArea_DEM',
    folder: EXPORT_FOLDER,
    scale: MAP_SCALE,
    region: roi,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF'
  };
  if (EXPORT_CRS) demParams.crs = EXPORT_CRS;

  Export.image.toDrive(demParams);
}

if (EXPORT_STUDYAREA_STATIC_TERRAIN_ZONES) {
  var staticMaps = terrainStack
    .addBands(zoneLow.rename('Zone_Low').toUint8())
    .addBands(zoneMid.rename('Zone_Mid').toUint8())
    .addBands(zoneHigh.rename('Zone_High').toUint8());

  var stParams = {
    image: staticMaps,
    description: 'Angren_Static_Terrain_Zones',
    folder: EXPORT_FOLDER,
    scale: MAP_SCALE,
    region: roi,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF'
  };
  if (EXPORT_CRS) stParams.crs = EXPORT_CRS;

  Export.image.toDrive(stParams);
}

// ----------------------------
// 2) FULL MAP EXPORTS: each Year × Season (NDVI, LST_C, ValidMask)
// ----------------------------
function makeComposite(year, seasonName) {
  var m = seasonMonths[seasonName];

  var col = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
    .merge(ee.ImageCollection('LANDSAT/LE07/C02/T1_L2'))
    .merge(ee.ImageCollection('LANDSAT/LC08/C02/T1_L2'))
    .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2'))
    .filterBounds(roi)
    .filter(ee.Filter.calendarRange(year, year, 'year'))
    .filter(ee.Filter.calendarRange(m.start, m.end, 'month'))
    .map(maskLandsatStrict)
    .map(applyScaleFactors)
    .map(unifyBands)
    .map(addIndices);

  return col.median().clip(roi);
}

function exportAllYearSeasonMapsFromFinalData() {
  var fc = finalData;

  if (MAPS_ONLY_WHEN_DATA_EXISTS) {
    fc = fc.filter(ee.Filter.gte('Image_Count', MAPS_MIN_SCENES));
  }
  if (MAPS_QC_GATE) {
    fc = fc.filter(ee.Filter.gte('Cov_Composite', MAPS_MIN_COV));
  }

  fc = fc.select(['Year', 'Season', 'Image_Count', 'Cov_Composite']);

  fc.toList(fc.size()).evaluate(function(list) {
    if (!list || list.length === 0) {
      print('No year-season combinations passed the export filters.');
      return;
    }

    print('Year-season exports to be created:', list.length);

    list.forEach(function(f) {
      var p = f.properties;
      var year = Number(p.Year);
      var seasonName = String(p.Season);

      if (!seasonMonths[seasonName]) return;

      var comp = makeComposite(year, seasonName);
      var tag = seasonName + '_' + year;

      if (EXPORT_NDVI) {
        var ndviParams = {
          image: comp.select('NDVI').toFloat(),
          description: 'Angren_NDVI_' + tag,
          folder: EXPORT_FOLDER,
          scale: MAP_SCALE,
          region: roi,
          maxPixels: 1e13,
          fileFormat: 'GeoTIFF'
        };
        if (EXPORT_CRS) ndviParams.crs = EXPORT_CRS;
        Export.image.toDrive(ndviParams);
      }

      if (EXPORT_LST) {
        var lstParams = {
          image: comp.select('LST_C').toFloat(),
          description: 'Angren_LST_' + tag,
          folder: EXPORT_FOLDER,
          scale: MAP_SCALE,
          region: roi,
          maxPixels: 1e13,
          fileFormat: 'GeoTIFF'
        };
        if (EXPORT_CRS) lstParams.crs = EXPORT_CRS;
        Export.image.toDrive(lstParams);
      }

      if (EXPORT_MASK) {
        var maskParams = {
          image: comp.select('Red').mask().rename('ValidMask').toUint8(),
          description: 'Angren_ValidMask_' + tag,
          folder: EXPORT_FOLDER,
          scale: MAP_SCALE,
          region: roi,
          maxPixels: 1e13,
          fileFormat: 'GeoTIFF'
        };
        if (EXPORT_CRS) maskParams.crs = EXPORT_CRS;
        Export.image.toDrive(maskParams);
      }
    });

    print('All export tasks created. Open the Tasks tab and run them.');
  });
}

// =====================================================================
// FULL EXPORT (PER YEAR × PER SEASON) + STUDY AREA PACK — FIXED CRS
// Key fix: DO NOT pass ee.Projection (refProj) into Export "crs".
// Use a CRS STRING (recommended) or omit "crs" entirely.
// =====================================================================

// ----------------------------
// 0) EXPORT SETTINGS (edit here)
// ----------------------------
var EXPORT_FOLDER = 'Angren_JournalPack';
var MAP_SCALE = 30;

// ✅ Recommended fixed CRS for Angren (~70E, 41N) = UTM 42N
// If you prefer to keep native image projection, set EXPORT_CRS = null
var EXPORT_CRS = 'EPSG:32642';   // or set to null to omit crs

var MAPS_ONLY_WHEN_DATA_EXISTS = true;

var MAPS_QC_GATE = true;
var MAPS_MIN_COV = 0.80;
var MAPS_MIN_SCENES = 1;

var EXPORT_NDVI = true;
var EXPORT_LST  = true;
var EXPORT_MASK = true;

var EXPORT_STUDYAREA_VECTOR = true;
var EXPORT_STUDYAREA_HILLSHADE = true;
var EXPORT_STUDYAREA_DEM = true;
var EXPORT_STUDYAREA_STATIC_TERRAIN_ZONES = true;

var seasonMonths = {
  'Spring': {start: 3, end: 5},
  'Summer': {start: 6, end: 8},
  'Autumn': {start: 9, end: 11}
};

// ----------------------------
// 1) STUDY AREA PACK (Figure 1 base)
// ----------------------------
if (EXPORT_STUDYAREA_VECTOR) {
  var roiFc = ee.FeatureCollection([ee.Feature(roi, {'name': 'ROI_Angren'})]);
  Export.table.toDrive({
    collection: roiFc,
    description: 'Angren_StudyArea_ROI',
    folder: EXPORT_FOLDER,
    fileFormat: 'SHP'
  });
}

if (EXPORT_STUDYAREA_HILLSHADE) {
  var hillshade = ee.Terrain.hillshade(elevation).rename('Hillshade');

  var hsParams = {
    image: hillshade.toUint8(),
    description: 'Angren_StudyArea_Hillshade',
    folder: EXPORT_FOLDER,
    scale: MAP_SCALE,
    region: roi,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF'
  };
  if (EXPORT_CRS) hsParams.crs = EXPORT_CRS;

  Export.image.toDrive(hsParams);
}

if (EXPORT_STUDYAREA_DEM) {
  var demParams = {
    image: elevation.rename('Elevation_m').toFloat(),
    description: 'Angren_StudyArea_DEM',
    folder: EXPORT_FOLDER,
    scale: MAP_SCALE,
    region: roi,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF'
  };
  if (EXPORT_CRS) demParams.crs = EXPORT_CRS;

  Export.image.toDrive(demParams);
}

if (EXPORT_STUDYAREA_STATIC_TERRAIN_ZONES) {
  var staticMaps = terrainStack
    .addBands(zoneLow.rename('Zone_Low').toUint8())
    .addBands(zoneMid.rename('Zone_Mid').toUint8())
    .addBands(zoneHigh.rename('Zone_High').toUint8());

  var stParams = {
    image: staticMaps,
    description: 'Angren_Static_Terrain_Zones',
    folder: EXPORT_FOLDER,
    scale: MAP_SCALE,
    region: roi,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF'
  };
  if (EXPORT_CRS) stParams.crs = EXPORT_CRS;

  Export.image.toDrive(stParams);
}

// ----------------------------
// 2) FULL MAP EXPORTS: each Year × Season (NDVI, LST_C, ValidMask)
// ----------------------------
function makeComposite(year, seasonName) {
  var m = seasonMonths[seasonName];

  var col = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
    .merge(ee.ImageCollection('LANDSAT/LE07/C02/T1_L2'))
    .merge(ee.ImageCollection('LANDSAT/LC08/C02/T1_L2'))
    .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2'))
    .filterBounds(roi)
    .filter(ee.Filter.calendarRange(year, year, 'year'))
    .filter(ee.Filter.calendarRange(m.start, m.end, 'month'))
    .map(maskLandsatStrict)
    .map(applyScaleFactors)
    .map(unifyBands)
    .map(addIndices);

  return col.median().clip(roi);
}

function exportAllYearSeasonMapsFromFinalData() {
  var fc = finalData;

  if (MAPS_ONLY_WHEN_DATA_EXISTS) {
    fc = fc.filter(ee.Filter.gte('Image_Count', MAPS_MIN_SCENES));
  }
  if (MAPS_QC_GATE) {
    fc = fc.filter(ee.Filter.gte('Cov_Composite', MAPS_MIN_COV));
  }

  fc = fc.select(['Year', 'Season', 'Image_Count', 'Cov_Composite']);

  fc.toList(fc.size()).evaluate(function(list) {
    if (!list || list.length === 0) {
      print('No year-season combinations passed the export filters.');
      return;
    }

    print('Year-season exports to be created:', list.length);

    list.forEach(function(f) {
      var p = f.properties;
      var year = Number(p.Year);
      var seasonName = String(p.Season);

      if (!seasonMonths[seasonName]) return;

      var comp = makeComposite(year, seasonName);
      var tag = seasonName + '_' + year;

      if (EXPORT_NDVI) {
        var ndviParams = {
          image: comp.select('NDVI').toFloat(),
          description: 'Angren_NDVI_' + tag,
          folder: EXPORT_FOLDER,
          scale: MAP_SCALE,
          region: roi,
          maxPixels: 1e13,
          fileFormat: 'GeoTIFF'
        };
        if (EXPORT_CRS) ndviParams.crs = EXPORT_CRS;
        Export.image.toDrive(ndviParams);
      }

      if (EXPORT_LST) {
        var lstParams = {
          image: comp.select('LST_C').toFloat(),
          description: 'Angren_LST_' + tag,
          folder: EXPORT_FOLDER,
          scale: MAP_SCALE,
          region: roi,
          maxPixels: 1e13,
          fileFormat: 'GeoTIFF'
        };
        if (EXPORT_CRS) lstParams.crs = EXPORT_CRS;
        Export.image.toDrive(lstParams);
      }

      if (EXPORT_MASK) {
        var maskParams = {
          image: comp.select('Red').mask().rename('ValidMask').toUint8(),
          description: 'Angren_ValidMask_' + tag,
          folder: EXPORT_FOLDER,
          scale: MAP_SCALE,
          region: roi,
          maxPixels: 1e13,
          fileFormat: 'GeoTIFF'
        };
        if (EXPORT_CRS) maskParams.crs = EXPORT_CRS;
        Export.image.toDrive(maskParams);
      }
    });

    print('All export tasks created. Open the Tasks tab and run them.');
  });
}

// =====================================================================
// FULL EXPORT (PER YEAR × PER SEASON) + STUDY AREA PACK — FIXED CRS
// Key fix: DO NOT pass ee.Projection (refProj) into Export "crs".
// Use a CRS STRING (recommended) or omit "crs" entirely.
// =====================================================================

// ----------------------------
// 0) EXPORT SETTINGS (edit here)
// ----------------------------
var EXPORT_FOLDER = 'Angren_JournalPack';
var MAP_SCALE = 30;

// ✅ Recommended fixed CRS for Angren (~70E, 41N) = UTM 42N
// If you prefer to keep native image projection, set EXPORT_CRS = null
var EXPORT_CRS = 'EPSG:32642';   // or set to null to omit crs

var MAPS_ONLY_WHEN_DATA_EXISTS = true;

var MAPS_QC_GATE = true;
var MAPS_MIN_COV = 0.80;
var MAPS_MIN_SCENES = 1;

var EXPORT_NDVI = true;
var EXPORT_LST  = true;
var EXPORT_MASK = true;

var EXPORT_STUDYAREA_VECTOR = true;
var EXPORT_STUDYAREA_HILLSHADE = true;
var EXPORT_STUDYAREA_DEM = true;
var EXPORT_STUDYAREA_STATIC_TERRAIN_ZONES = true;

var seasonMonths = {
  'Spring': {start: 3, end: 5},
  'Summer': {start: 6, end: 8},
  'Autumn': {start: 9, end: 11}
};

// ----------------------------
// 1) STUDY AREA PACK (Figure 1 base)
// ----------------------------
if (EXPORT_STUDYAREA_VECTOR) {
  var roiFc = ee.FeatureCollection([ee.Feature(roi, {'name': 'ROI_Angren'})]);
  Export.table.toDrive({
    collection: roiFc,
    description: 'Angren_StudyArea_ROI',
    folder: EXPORT_FOLDER,
    fileFormat: 'SHP'
  });
}

if (EXPORT_STUDYAREA_HILLSHADE) {
  var hillshade = ee.Terrain.hillshade(elevation).rename('Hillshade');

  var hsParams = {
    image: hillshade.toUint8(),
    description: 'Angren_StudyArea_Hillshade',
    folder: EXPORT_FOLDER,
    scale: MAP_SCALE,
    region: roi,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF'
  };
  if (EXPORT_CRS) hsParams.crs = EXPORT_CRS;

  Export.image.toDrive(hsParams);
}

if (EXPORT_STUDYAREA_DEM) {
  var demParams = {
    image: elevation.rename('Elevation_m').toFloat(),
    description: 'Angren_StudyArea_DEM',
    folder: EXPORT_FOLDER,
    scale: MAP_SCALE,
    region: roi,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF'
  };
  if (EXPORT_CRS) demParams.crs = EXPORT_CRS;

  Export.image.toDrive(demParams);
}

if (EXPORT_STUDYAREA_STATIC_TERRAIN_ZONES) {
  var staticMaps = terrainStack
    .addBands(zoneLow.rename('Zone_Low').toUint8())
    .addBands(zoneMid.rename('Zone_Mid').toUint8())
    .addBands(zoneHigh.rename('Zone_High').toUint8());

  var stParams = {
    image: staticMaps,
    description: 'Angren_Static_Terrain_Zones',
    folder: EXPORT_FOLDER,
    scale: MAP_SCALE,
    region: roi,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF'
  };
  if (EXPORT_CRS) stParams.crs = EXPORT_CRS;

  Export.image.toDrive(stParams);
}

// ----------------------------
// 2) FULL MAP EXPORTS: each Year × Season (NDVI, LST_C, ValidMask)
// ----------------------------
function makeComposite(year, seasonName) {
  var m = seasonMonths[seasonName];

  var col = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
    .merge(ee.ImageCollection('LANDSAT/LE07/C02/T1_L2'))
    .merge(ee.ImageCollection('LANDSAT/LC08/C02/T1_L2'))
    .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2'))
    .filterBounds(roi)
    .filter(ee.Filter.calendarRange(year, year, 'year'))
    .filter(ee.Filter.calendarRange(m.start, m.end, 'month'))
    .map(maskLandsatStrict)
    .map(applyScaleFactors)
    .map(unifyBands)
    .map(addIndices);

  return col.median().clip(roi);
}

function exportAllYearSeasonMapsFromFinalData() {
  var fc = finalData;

  if (MAPS_ONLY_WHEN_DATA_EXISTS) {
    fc = fc.filter(ee.Filter.gte('Image_Count', MAPS_MIN_SCENES));
  }
  if (MAPS_QC_GATE) {
    fc = fc.filter(ee.Filter.gte('Cov_Composite', MAPS_MIN_COV));
  }

  fc = fc.select(['Year', 'Season', 'Image_Count', 'Cov_Composite']);

  fc.toList(fc.size()).evaluate(function(list) {
    if (!list || list.length === 0) {
      print('No year-season combinations passed the export filters.');
      return;
    }

    print('Year-season exports to be created:', list.length);

    list.forEach(function(f) {
      var p = f.properties;
      var year = Number(p.Year);
      var seasonName = String(p.Season);

      if (!seasonMonths[seasonName]) return;

      var comp = makeComposite(year, seasonName);
      var tag = seasonName + '_' + year;

      if (EXPORT_NDVI) {
        var ndviParams = {
          image: comp.select('NDVI').toFloat(),
          description: 'Angren_NDVI_' + tag,
          folder: EXPORT_FOLDER,
          scale: MAP_SCALE,
          region: roi,
          maxPixels: 1e13,
          fileFormat: 'GeoTIFF'
        };
        if (EXPORT_CRS) ndviParams.crs = EXPORT_CRS;
        Export.image.toDrive(ndviParams);
      }

      if (EXPORT_LST) {
        var lstParams = {
          image: comp.select('LST_C').toFloat(),
          description: 'Angren_LST_' + tag,
          folder: EXPORT_FOLDER,
          scale: MAP_SCALE,
          region: roi,
          maxPixels: 1e13,
          fileFormat: 'GeoTIFF'
        };
        if (EXPORT_CRS) lstParams.crs = EXPORT_CRS;
        Export.image.toDrive(lstParams);
      }

      if (EXPORT_MASK) {
        var maskParams = {
          image: comp.select('Red').mask().rename('ValidMask').toUint8(),
          description: 'Angren_ValidMask_' + tag,
          folder: EXPORT_FOLDER,
          scale: MAP_SCALE,
          region: roi,
          maxPixels: 1e13,
          fileFormat: 'GeoTIFF'
        };
        if (EXPORT_CRS) maskParams.crs = EXPORT_CRS;
        Export.image.toDrive(maskParams);
      }
    });

    print('All export tasks created. Open the Tasks tab and run them.');
  });
}

// =====================================================================
// FULL EXPORT (PER YEAR × PER SEASON) + STUDY AREA PACK — FIXED CRS
// Key fix: DO NOT pass ee.Projection (refProj) into Export "crs".
// Use a CRS STRING (recommended) or omit "crs" entirely.
// =====================================================================

// ----------------------------
// 0) EXPORT SETTINGS (edit here)
// ----------------------------
var EXPORT_FOLDER = 'Angren_JournalPack';
var MAP_SCALE = 30;

// ✅ Recommended fixed CRS for Angren (~70E, 41N) = UTM 42N
// If you prefer to keep native image projection, set EXPORT_CRS = null
var EXPORT_CRS = 'EPSG:32642';   // or set to null to omit crs

var MAPS_ONLY_WHEN_DATA_EXISTS = true;

var MAPS_QC_GATE = true;
var MAPS_MIN_COV = 0.80;
var MAPS_MIN_SCENES = 1;

var EXPORT_NDVI = true;
var EXPORT_LST  = true;
var EXPORT_MASK = true;

var EXPORT_STUDYAREA_VECTOR = true;
var EXPORT_STUDYAREA_HILLSHADE = true;
var EXPORT_STUDYAREA_DEM = true;
var EXPORT_STUDYAREA_STATIC_TERRAIN_ZONES = true;

var seasonMonths = {
  'Spring': {start: 3, end: 5},
  'Summer': {start: 6, end: 8},
  'Autumn': {start: 9, end: 11}
};

// ----------------------------
// 1) STUDY AREA PACK (Figure 1 base)
// ----------------------------
if (EXPORT_STUDYAREA_VECTOR) {
  var roiFc = ee.FeatureCollection([ee.Feature(roi, {'name': 'ROI_Angren'})]);
  Export.table.toDrive({
    collection: roiFc,
    description: 'Angren_StudyArea_ROI',
    folder: EXPORT_FOLDER,
    fileFormat: 'SHP'
  });
}

if (EXPORT_STUDYAREA_HILLSHADE) {
  var hillshade = ee.Terrain.hillshade(elevation).rename('Hillshade');

  var hsParams = {
    image: hillshade.toUint8(),
    description: 'Angren_StudyArea_Hillshade',
    folder: EXPORT_FOLDER,
    scale: MAP_SCALE,
    region: roi,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF'
  };
  if (EXPORT_CRS) hsParams.crs = EXPORT_CRS;

  Export.image.toDrive(hsParams);
}

if (EXPORT_STUDYAREA_DEM) {
  var demParams = {
    image: elevation.rename('Elevation_m').toFloat(),
    description: 'Angren_StudyArea_DEM',
    folder: EXPORT_FOLDER,
    scale: MAP_SCALE,
    region: roi,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF'
  };
  if (EXPORT_CRS) demParams.crs = EXPORT_CRS;

  Export.image.toDrive(demParams);
}

if (EXPORT_STUDYAREA_STATIC_TERRAIN_ZONES) {
  var staticMaps = terrainStack
    .addBands(zoneLow.rename('Zone_Low').toUint8())
    .addBands(zoneMid.rename('Zone_Mid').toUint8())
    .addBands(zoneHigh.rename('Zone_High').toUint8());

  var stParams = {
    image: staticMaps,
    description: 'Angren_Static_Terrain_Zones',
    folder: EXPORT_FOLDER,
    scale: MAP_SCALE,
    region: roi,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF'
  };
  if (EXPORT_CRS) stParams.crs = EXPORT_CRS;

  Export.image.toDrive(stParams);
}

// ----------------------------
// 2) FULL MAP EXPORTS: each Year × Season (NDVI, LST_C, ValidMask)
// ----------------------------
function makeComposite(year, seasonName) {
  var m = seasonMonths[seasonName];

  var col = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
    .merge(ee.ImageCollection('LANDSAT/LE07/C02/T1_L2'))
    .merge(ee.ImageCollection('LANDSAT/LC08/C02/T1_L2'))
    .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2'))
    .filterBounds(roi)
    .filter(ee.Filter.calendarRange(year, year, 'year'))
    .filter(ee.Filter.calendarRange(m.start, m.end, 'month'))
    .map(maskLandsatStrict)
    .map(applyScaleFactors)
    .map(unifyBands)
    .map(addIndices);

  return col.median().clip(roi);
}

function exportAllYearSeasonMapsFromFinalData() {
  var fc = finalData;

  if (MAPS_ONLY_WHEN_DATA_EXISTS) {
    fc = fc.filter(ee.Filter.gte('Image_Count', MAPS_MIN_SCENES));
  }
  if (MAPS_QC_GATE) {
    fc = fc.filter(ee.Filter.gte('Cov_Composite', MAPS_MIN_COV));
  }

  fc = fc.select(['Year', 'Season', 'Image_Count', 'Cov_Composite']);

  fc.toList(fc.size()).evaluate(function(list) {
    if (!list || list.length === 0) {
      print('No year-season combinations passed the export filters.');
      return;
    }

    print('Year-season exports to be created:', list.length);

    list.forEach(function(f) {
      var p = f.properties;
      var year = Number(p.Year);
      var seasonName = String(p.Season);

      if (!seasonMonths[seasonName]) return;

      var comp = makeComposite(year, seasonName);
      var tag = seasonName + '_' + year;

      if (EXPORT_NDVI) {
        var ndviParams = {
          image: comp.select('NDVI').toFloat(),
          description: 'Angren_NDVI_' + tag,
          folder: EXPORT_FOLDER,
          scale: MAP_SCALE,
          region: roi,
          maxPixels: 1e13,
          fileFormat: 'GeoTIFF'
        };
        if (EXPORT_CRS) ndviParams.crs = EXPORT_CRS;
        Export.image.toDrive(ndviParams);
      }

      if (EXPORT_LST) {
        var lstParams = {
          image: comp.select('LST_C').toFloat(),
          description: 'Angren_LST_' + tag,
          folder: EXPORT_FOLDER,
          scale: MAP_SCALE,
          region: roi,
          maxPixels: 1e13,
          fileFormat: 'GeoTIFF'
        };
        if (EXPORT_CRS) lstParams.crs = EXPORT_CRS;
        Export.image.toDrive(lstParams);
      }

      if (EXPORT_MASK) {
        var maskParams = {
          image: comp.select('Red').mask().rename('ValidMask').toUint8(),
          description: 'Angren_ValidMask_' + tag,
          folder: EXPORT_FOLDER,
          scale: MAP_SCALE,
          region: roi,
          maxPixels: 1e13,
          fileFormat: 'GeoTIFF'
        };
        if (EXPORT_CRS) maskParams.crs = EXPORT_CRS;
        Export.image.toDrive(maskParams);
      }
    });

    print('All export tasks created. Open the Tasks tab and run them.');
  });
}



// =====================================================================
// FULL EXPORT (PER YEAR × PER SEASON) + STUDY AREA PACK
// Paste this AFTER you define: finalData, refProj, roi, elevation, terrainStack,
// and functions: maskLandsatStrict, applyScaleFactors, unifyBands, addIndices
// =====================================================================

// ----------------------------
// 0) EXPORT SETTINGS (edit here)
// ----------------------------
var EXPORT_FOLDER = 'Angren_JournalPack';   // Drive folder (will be created if absent)
var MAP_SCALE = 30;

// Export only valid year-season combos from finalData (recommended)
var MAPS_ONLY_WHEN_DATA_EXISTS = true;

// Optional stricter gate for map exports (recommended for journal figures)
var MAPS_QC_GATE = true;
var MAPS_MIN_COV = 0.80;     // export maps only if Cov_Composite >= this
var MAPS_MIN_SCENES = 1;     // export maps only if Image_Count >= this

// What to export per year-season
var EXPORT_NDVI = true;
var EXPORT_LST  = true;
var EXPORT_MASK = true;

// Study area exports
var EXPORT_STUDYAREA_VECTOR = true;
var EXPORT_STUDYAREA_HILLSHADE = true;
var EXPORT_STUDYAREA_DEM = true;
var EXPORT_STUDYAREA_STATIC_TERRAIN_ZONES = true;

// Season month lookup (must match your seasons)
var seasonMonths = {
  'Spring': {start: 3, end: 5},
  'Summer': {start: 6, end: 8},
  'Autumn': {start: 9, end: 11}
};

// ----------------------------
// 1) STUDY AREA PACK (Figure 1 base)
// ----------------------------
if (EXPORT_STUDYAREA_VECTOR) {
  var roiFc = ee.FeatureCollection([ee.Feature(roi, {'name': 'ROI_Angren'})]);
  Export.table.toDrive({
    collection: roiFc,
    description: 'Angren_StudyArea_ROI',
    folder: EXPORT_FOLDER,
    fileFormat: 'SHP'
  });
}

if (EXPORT_STUDYAREA_HILLSHADE) {
  var hillshade = ee.Terrain.hillshade(elevation).rename('Hillshade');
  Export.image.toDrive({
    image: hillshade.toUint8(),
    description: 'Angren_StudyArea_Hillshade',
    folder: EXPORT_FOLDER,
    scale: MAP_SCALE,
    region: roi,
    crs: refProj,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF'
  });
}

if (EXPORT_STUDYAREA_DEM) {
  Export.image.toDrive({
    image: elevation.rename('Elevation_m').toFloat(),
    description: 'Angren_StudyArea_DEM',
    folder: EXPORT_FOLDER,
    scale: MAP_SCALE,
    region: roi,
    crs: refProj,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF'
  });
}

if (EXPORT_STUDYAREA_STATIC_TERRAIN_ZONES) {
  // Requires zoneLow/zoneMid/zoneHigh already defined in your script
  var staticMaps = terrainStack
    .addBands(zoneLow.rename('Zone_Low').toUint8())
    .addBands(zoneMid.rename('Zone_Mid').toUint8())
    .addBands(zoneHigh.rename('Zone_High').toUint8());

  Export.image.toDrive({
    image: staticMaps,
    description: 'Angren_Static_Terrain_Zones',
    folder: EXPORT_FOLDER,
    scale: MAP_SCALE,
    region: roi,
    crs: refProj,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF'
  });
}

// ----------------------------
// 2) FULL MAP EXPORTS: each Year × Season (NDVI, LST_C, ValidMask)
// ----------------------------

// Helper: build composite for a given year + seasonName
function makeComposite(year, seasonName) {
  var m = seasonMonths[seasonName];

  var col = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
    .merge(ee.ImageCollection('LANDSAT/LE07/C02/T1_L2'))
    .merge(ee.ImageCollection('LANDSAT/LC08/C02/T1_L2'))
    .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2'))
    .filterBounds(roi)
    .filter(ee.Filter.calendarRange(year, year, 'year'))
    .filter(ee.Filter.calendarRange(m.start, m.end, 'month'))
    .map(maskLandsatStrict)
    .map(applyScaleFactors)
    .map(unifyBands)
    .map(addIndices);

  return col.median().clip(roi);
}

// Create exports only for combos that exist in finalData (no empty years)
function exportAllYearSeasonMapsFromFinalData() {
  var fc = finalData;

  if (MAPS_ONLY_WHEN_DATA_EXISTS) {
    fc = fc.filter(ee.Filter.gte('Image_Count', MAPS_MIN_SCENES));
  }
  if (MAPS_QC_GATE) {
    fc = fc.filter(ee.Filter.gte('Cov_Composite', MAPS_MIN_COV));
  }

  // Keep only what we need client-side
  fc = fc.select(['Year', 'Season', 'Image_Count', 'Cov_Composite']);

  // Bring list of features to client side, then create Export tasks
  fc.toList(fc.size()).evaluate(function(list) {
    if (!list || list.length === 0) {
      print('No year-season combinations passed the export filters.');
      return;
    }

    print('Year-season exports to be created:', list.length);

    list.forEach(function(f) {
      var p = f.properties;
      var year = Number(p.Year);
      var seasonName = String(p.Season);

      // Safety: only export for known seasons
      if (!seasonMonths[seasonName]) return;

      var comp = makeComposite(year, seasonName);

      var tag = seasonName + '_' + year;

      if (EXPORT_NDVI) {
        Export.image.toDrive({
          image: comp.select('NDVI').toFloat(),
          description: 'Angren_NDVI_' + tag,
          folder: EXPORT_FOLDER,
          scale: MAP_SCALE,
          region: roi,
          crs: refProj,
          maxPixels: 1e13,
          fileFormat: 'GeoTIFF'
        });
      }

      if (EXPORT_LST) {
        Export.image.toDrive({
          image: comp.select('LST_C').toFloat(),
          description: 'Angren_LST_' + tag,
          folder: EXPORT_FOLDER,
          scale: MAP_SCALE,
          region: roi,
          crs: refProj,
          maxPixels: 1e13,
          fileFormat: 'GeoTIFF'
        });
      }

      if (EXPORT_MASK) {
        Export.image.toDrive({
          image: comp.select('Red').mask().rename('ValidMask').toUint8(),
          description: 'Angren_ValidMask_' + tag,
          folder: EXPORT_FOLDER,
          scale: MAP_SCALE,
          region: roi,
          crs: refProj,
          maxPixels: 1e13,
          fileFormat: 'GeoTIFF'
        });
      }
    });

    print('All export tasks created. Open the Tasks tab and run them.');
  });
}


// =====================================================================
// FULL EXPORT (PER YEAR × PER SEASON) + STUDY AREA PACK — FIXED CRS
// Key fix: DO NOT pass ee.Projection (refProj) into Export "crs".
// Use a CRS STRING (recommended) or omit "crs" entirely.
// =====================================================================

// ----------------------------
// 0) EXPORT SETTINGS (edit here)
// ----------------------------
var EXPORT_FOLDER = 'Angren_JournalPack';
var MAP_SCALE = 30;

// ✅ Recommended fixed CRS for Angren (~70E, 41N) = UTM 42N
// If you prefer to keep native image projection, set EXPORT_CRS = null
var EXPORT_CRS = 'EPSG:32642';   // or set to null to omit crs

var MAPS_ONLY_WHEN_DATA_EXISTS = true;

var MAPS_QC_GATE = true;
var MAPS_MIN_COV = 0.80;
var MAPS_MIN_SCENES = 1;

var EXPORT_NDVI = true;
var EXPORT_LST  = true;
var EXPORT_MASK = true;

var EXPORT_STUDYAREA_VECTOR = true;
var EXPORT_STUDYAREA_HILLSHADE = true;
var EXPORT_STUDYAREA_DEM = true;
var EXPORT_STUDYAREA_STATIC_TERRAIN_ZONES = true;

var seasonMonths = {
  'Spring': {start: 3, end: 5},
  'Summer': {start: 6, end: 8},
  'Autumn': {start: 9, end: 11}
};

// ----------------------------
// 1) STUDY AREA PACK (Figure 1 base)
// ----------------------------
if (EXPORT_STUDYAREA_VECTOR) {
  var roiFc = ee.FeatureCollection([ee.Feature(roi, {'name': 'ROI_Angren'})]);
  Export.table.toDrive({
    collection: roiFc,
    description: 'Angren_StudyArea_ROI',
    folder: EXPORT_FOLDER,
    fileFormat: 'SHP'
  });
}

if (EXPORT_STUDYAREA_HILLSHADE) {
  var hillshade = ee.Terrain.hillshade(elevation).rename('Hillshade');

  var hsParams = {
    image: hillshade.toUint8(),
    description: 'Angren_StudyArea_Hillshade',
    folder: EXPORT_FOLDER,
    scale: MAP_SCALE,
    region: roi,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF'
  };
  if (EXPORT_CRS) hsParams.crs = EXPORT_CRS;

  Export.image.toDrive(hsParams);
}

if (EXPORT_STUDYAREA_DEM) {
  var demParams = {
    image: elevation.rename('Elevation_m').toFloat(),
    description: 'Angren_StudyArea_DEM',
    folder: EXPORT_FOLDER,
    scale: MAP_SCALE,
    region: roi,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF'
  };
  if (EXPORT_CRS) demParams.crs = EXPORT_CRS;

  Export.image.toDrive(demParams);
}

if (EXPORT_STUDYAREA_STATIC_TERRAIN_ZONES) {
  var staticMaps = terrainStack
    .addBands(zoneLow.rename('Zone_Low').toUint8())
    .addBands(zoneMid.rename('Zone_Mid').toUint8())
    .addBands(zoneHigh.rename('Zone_High').toUint8());

  var stParams = {
    image: staticMaps,
    description: 'Angren_Static_Terrain_Zones',
    folder: EXPORT_FOLDER,
    scale: MAP_SCALE,
    region: roi,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF'
  };
  if (EXPORT_CRS) stParams.crs = EXPORT_CRS;

  Export.image.toDrive(stParams);
}

// ----------------------------
// 2) FULL MAP EXPORTS: each Year × Season (NDVI, LST_C, ValidMask)
// ----------------------------
function makeComposite(year, seasonName) {
  var m = seasonMonths[seasonName];

  var col = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
    .merge(ee.ImageCollection('LANDSAT/LE07/C02/T1_L2'))
    .merge(ee.ImageCollection('LANDSAT/LC08/C02/T1_L2'))
    .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2'))
    .filterBounds(roi)
    .filter(ee.Filter.calendarRange(year, year, 'year'))
    .filter(ee.Filter.calendarRange(m.start, m.end, 'month'))
    .map(maskLandsatStrict)
    .map(applyScaleFactors)
    .map(unifyBands)
    .map(addIndices);

  return col.median().clip(roi);
}

function exportAllYearSeasonMapsFromFinalData() {
  var fc = finalData;

  if (MAPS_ONLY_WHEN_DATA_EXISTS) {
    fc = fc.filter(ee.Filter.gte('Image_Count', MAPS_MIN_SCENES));
  }
  if (MAPS_QC_GATE) {
    fc = fc.filter(ee.Filter.gte('Cov_Composite', MAPS_MIN_COV));
  }

  fc = fc.select(['Year', 'Season', 'Image_Count', 'Cov_Composite']);

  fc.toList(fc.size()).evaluate(function(list) {
    if (!list || list.length === 0) {
      print('No year-season combinations passed the export filters.');
      return;
    }

    print('Year-season exports to be created:', list.length);

    list.forEach(function(f) {
      var p = f.properties;
      var year = Number(p.Year);
      var seasonName = String(p.Season);

      if (!seasonMonths[seasonName]) return;

      var comp = makeComposite(year, seasonName);
      var tag = seasonName + '_' + year;

      if (EXPORT_NDVI) {
        var ndviParams = {
          image: comp.select('NDVI').toFloat(),
          description: 'Angren_NDVI_' + tag,
          folder: EXPORT_FOLDER,
          scale: MAP_SCALE,
          region: roi,
          maxPixels: 1e13,
          fileFormat: 'GeoTIFF'
        };
        if (EXPORT_CRS) ndviParams.crs = EXPORT_CRS;
        Export.image.toDrive(ndviParams);
      }

      if (EXPORT_LST) {
        var lstParams = {
          image: comp.select('LST_C').toFloat(),
          description: 'Angren_LST_' + tag,
          folder: EXPORT_FOLDER,
          scale: MAP_SCALE,
          region: roi,
          maxPixels: 1e13,
          fileFormat: 'GeoTIFF'
        };
        if (EXPORT_CRS) lstParams.crs = EXPORT_CRS;
        Export.image.toDrive(lstParams);
      }

      if (EXPORT_MASK) {
        var maskParams = {
          image: comp.select('Red').mask().rename('ValidMask').toUint8(),
          description: 'Angren_ValidMask_' + tag,
          folder: EXPORT_FOLDER,
          scale: MAP_SCALE,
          region: roi,
          maxPixels: 1e13,
          fileFormat: 'GeoTIFF'
        };
        if (EXPORT_CRS) maskParams.crs = EXPORT_CRS;
        Export.image.toDrive(maskParams);
      }
    });

    print('All export tasks created. Open the Tasks tab and run them.');
  });
}

// =====================================================================
// FULL EXPORT (PER YEAR × PER SEASON) + STUDY AREA PACK — FIXED CRS
// Key fix: DO NOT pass ee.Projection (refProj) into Export "crs".
// Use a CRS STRING (recommended) or omit "crs" entirely.
// =====================================================================

// ----------------------------
// 0) EXPORT SETTINGS (edit here)
// ----------------------------
var EXPORT_FOLDER = 'Angren_JournalPack';
var MAP_SCALE = 30;

// ✅ Recommended fixed CRS for Angren (~70E, 41N) = UTM 42N
// If you prefer to keep native image projection, set EXPORT_CRS = null
var EXPORT_CRS = 'EPSG:32642';   // or set to null to omit crs

var MAPS_ONLY_WHEN_DATA_EXISTS = true;

var MAPS_QC_GATE = true;
var MAPS_MIN_COV = 0.80;
var MAPS_MIN_SCENES = 1;

var EXPORT_NDVI = true;
var EXPORT_LST  = true;
var EXPORT_MASK = true;

var EXPORT_STUDYAREA_VECTOR = true;
var EXPORT_STUDYAREA_HILLSHADE = true;
var EXPORT_STUDYAREA_DEM = true;
var EXPORT_STUDYAREA_STATIC_TERRAIN_ZONES = true;

var seasonMonths = {
  'Spring': {start: 3, end: 5},
  'Summer': {start: 6, end: 8},
  'Autumn': {start: 9, end: 11}
};

// ----------------------------
// 1) STUDY AREA PACK (Figure 1 base)
// ----------------------------
if (EXPORT_STUDYAREA_VECTOR) {
  var roiFc = ee.FeatureCollection([ee.Feature(roi, {'name': 'ROI_Angren'})]);
  Export.table.toDrive({
    collection: roiFc,
    description: 'Angren_StudyArea_ROI',
    folder: EXPORT_FOLDER,
    fileFormat: 'SHP'
  });
}

if (EXPORT_STUDYAREA_HILLSHADE) {
  var hillshade = ee.Terrain.hillshade(elevation).rename('Hillshade');

  var hsParams = {
    image: hillshade.toUint8(),
    description: 'Angren_StudyArea_Hillshade',
    folder: EXPORT_FOLDER,
    scale: MAP_SCALE,
    region: roi,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF'
  };
  if (EXPORT_CRS) hsParams.crs = EXPORT_CRS;

  Export.image.toDrive(hsParams);
}

if (EXPORT_STUDYAREA_DEM) {
  var demParams = {
    image: elevation.rename('Elevation_m').toFloat(),
    description: 'Angren_StudyArea_DEM',
    folder: EXPORT_FOLDER,
    scale: MAP_SCALE,
    region: roi,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF'
  };
  if (EXPORT_CRS) demParams.crs = EXPORT_CRS;

  Export.image.toDrive(demParams);
}

if (EXPORT_STUDYAREA_STATIC_TERRAIN_ZONES) {
  var staticMaps = terrainStack
    .addBands(zoneLow.rename('Zone_Low').toUint8())
    .addBands(zoneMid.rename('Zone_Mid').toUint8())
    .addBands(zoneHigh.rename('Zone_High').toUint8());

  var stParams = {
    image: staticMaps,
    description: 'Angren_Static_Terrain_Zones',
    folder: EXPORT_FOLDER,
    scale: MAP_SCALE,
    region: roi,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF'
  };
  if (EXPORT_CRS) stParams.crs = EXPORT_CRS;

  Export.image.toDrive(stParams);
}

// ----------------------------
// 2) FULL MAP EXPORTS: each Year × Season (NDVI, LST_C, ValidMask)
// ----------------------------
function makeComposite(year, seasonName) {
  var m = seasonMonths[seasonName];

  var col = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
    .merge(ee.ImageCollection('LANDSAT/LE07/C02/T1_L2'))
    .merge(ee.ImageCollection('LANDSAT/LC08/C02/T1_L2'))
    .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2'))
    .filterBounds(roi)
    .filter(ee.Filter.calendarRange(year, year, 'year'))
    .filter(ee.Filter.calendarRange(m.start, m.end, 'month'))
    .map(maskLandsatStrict)
    .map(applyScaleFactors)
    .map(unifyBands)
    .map(addIndices);

  return col.median().clip(roi);
}

function exportAllYearSeasonMapsFromFinalData() {
  var fc = finalData;

  if (MAPS_ONLY_WHEN_DATA_EXISTS) {
    fc = fc.filter(ee.Filter.gte('Image_Count', MAPS_MIN_SCENES));
  }
  if (MAPS_QC_GATE) {
    fc = fc.filter(ee.Filter.gte('Cov_Composite', MAPS_MIN_COV));
  }

  fc = fc.select(['Year', 'Season', 'Image_Count', 'Cov_Composite']);

  fc.toList(fc.size()).evaluate(function(list) {
    if (!list || list.length === 0) {
      print('No year-season combinations passed the export filters.');
      return;
    }

    print('Year-season exports to be created:', list.length);

    list.forEach(function(f) {
      var p = f.properties;
      var year = Number(p.Year);
      var seasonName = String(p.Season);

      if (!seasonMonths[seasonName]) return;

      var comp = makeComposite(year, seasonName);
      var tag = seasonName + '_' + year;

      if (EXPORT_NDVI) {
        var ndviParams = {
          image: comp.select('NDVI').toFloat(),
          description: 'Angren_NDVI_' + tag,
          folder: EXPORT_FOLDER,
          scale: MAP_SCALE,
          region: roi,
          maxPixels: 1e13,
          fileFormat: 'GeoTIFF'
        };
        if (EXPORT_CRS) ndviParams.crs = EXPORT_CRS;
        Export.image.toDrive(ndviParams);
      }

      if (EXPORT_LST) {
        var lstParams = {
          image: comp.select('LST_C').toFloat(),
          description: 'Angren_LST_' + tag,
          folder: EXPORT_FOLDER,
          scale: MAP_SCALE,
          region: roi,
          maxPixels: 1e13,
          fileFormat: 'GeoTIFF'
        };
        if (EXPORT_CRS) lstParams.crs = EXPORT_CRS;
        Export.image.toDrive(lstParams);
      }

      if (EXPORT_MASK) {
        var maskParams = {
          image: comp.select('Red').mask().rename('ValidMask').toUint8(),
          description: 'Angren_ValidMask_' + tag,
          folder: EXPORT_FOLDER,
          scale: MAP_SCALE,
          region: roi,
          maxPixels: 1e13,
          fileFormat: 'GeoTIFF'
        };
        if (EXPORT_CRS) maskParams.crs = EXPORT_CRS;
        Export.image.toDrive(maskParams);
      }
    });

    print('All export tasks created. Open the Tasks tab and run them.');
  });
}

// =====================================================================
// FULL EXPORT (PER YEAR × PER SEASON) + STUDY AREA PACK — FIXED CRS
// Key fix: DO NOT pass ee.Projection (refProj) into Export "crs".
// Use a CRS STRING (recommended) or omit "crs" entirely.
// =====================================================================

// ----------------------------
// 0) EXPORT SETTINGS (edit here)
// ----------------------------
var EXPORT_FOLDER = 'Angren_JournalPack';
var MAP_SCALE = 30;

// ✅ Recommended fixed CRS for Angren (~70E, 41N) = UTM 42N
// If you prefer to keep native image projection, set EXPORT_CRS = null
var EXPORT_CRS = 'EPSG:32642';   // or set to null to omit crs

var MAPS_ONLY_WHEN_DATA_EXISTS = true;

var MAPS_QC_GATE = true;
var MAPS_MIN_COV = 0.80;
var MAPS_MIN_SCENES = 1;

var EXPORT_NDVI = true;
var EXPORT_LST  = true;
var EXPORT_MASK = true;

var EXPORT_STUDYAREA_VECTOR = true;
var EXPORT_STUDYAREA_HILLSHADE = true;
var EXPORT_STUDYAREA_DEM = true;
var EXPORT_STUDYAREA_STATIC_TERRAIN_ZONES = true;

var seasonMonths = {
  'Spring': {start: 3, end: 5},
  'Summer': {start: 6, end: 8},
  'Autumn': {start: 9, end: 11}
};

// ----------------------------
// 1) STUDY AREA PACK (Figure 1 base)
// ----------------------------
if (EXPORT_STUDYAREA_VECTOR) {
  var roiFc = ee.FeatureCollection([ee.Feature(roi, {'name': 'ROI_Angren'})]);
  Export.table.toDrive({
    collection: roiFc,
    description: 'Angren_StudyArea_ROI',
    folder: EXPORT_FOLDER,
    fileFormat: 'SHP'
  });
}

if (EXPORT_STUDYAREA_HILLSHADE) {
  var hillshade = ee.Terrain.hillshade(elevation).rename('Hillshade');

  var hsParams = {
    image: hillshade.toUint8(),
    description: 'Angren_StudyArea_Hillshade',
    folder: EXPORT_FOLDER,
    scale: MAP_SCALE,
    region: roi,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF'
  };
  if (EXPORT_CRS) hsParams.crs = EXPORT_CRS;

  Export.image.toDrive(hsParams);
}

if (EXPORT_STUDYAREA_DEM) {
  var demParams = {
    image: elevation.rename('Elevation_m').toFloat(),
    description: 'Angren_StudyArea_DEM',
    folder: EXPORT_FOLDER,
    scale: MAP_SCALE,
    region: roi,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF'
  };
  if (EXPORT_CRS) demParams.crs = EXPORT_CRS;

  Export.image.toDrive(demParams);
}

if (EXPORT_STUDYAREA_STATIC_TERRAIN_ZONES) {
  var staticMaps = terrainStack
    .addBands(zoneLow.rename('Zone_Low').toUint8())
    .addBands(zoneMid.rename('Zone_Mid').toUint8())
    .addBands(zoneHigh.rename('Zone_High').toUint8());

  var stParams = {
    image: staticMaps,
    description: 'Angren_Static_Terrain_Zones',
    folder: EXPORT_FOLDER,
    scale: MAP_SCALE,
    region: roi,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF'
  };
  if (EXPORT_CRS) stParams.crs = EXPORT_CRS;

  Export.image.toDrive(stParams);
}

// ----------------------------
// 2) FULL MAP EXPORTS: each Year × Season (NDVI, LST_C, ValidMask)
// ----------------------------
function makeComposite(year, seasonName) {
  var m = seasonMonths[seasonName];

  var col = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
    .merge(ee.ImageCollection('LANDSAT/LE07/C02/T1_L2'))
    .merge(ee.ImageCollection('LANDSAT/LC08/C02/T1_L2'))
    .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2'))
    .filterBounds(roi)
    .filter(ee.Filter.calendarRange(year, year, 'year'))
    .filter(ee.Filter.calendarRange(m.start, m.end, 'month'))
    .map(maskLandsatStrict)
    .map(applyScaleFactors)
    .map(unifyBands)
    .map(addIndices);

  return col.median().clip(roi);
}

function exportAllYearSeasonMapsFromFinalData() {
  var fc = finalData;

  if (MAPS_ONLY_WHEN_DATA_EXISTS) {
    fc = fc.filter(ee.Filter.gte('Image_Count', MAPS_MIN_SCENES));
  }
  if (MAPS_QC_GATE) {
    fc = fc.filter(ee.Filter.gte('Cov_Composite', MAPS_MIN_COV));
  }

  fc = fc.select(['Year', 'Season', 'Image_Count', 'Cov_Composite']);

  fc.toList(fc.size()).evaluate(function(list) {
    if (!list || list.length === 0) {
      print('No year-season combinations passed the export filters.');
      return;
    }

    print('Year-season exports to be created:', list.length);

    list.forEach(function(f) {
      var p = f.properties;
      var year = Number(p.Year);
      var seasonName = String(p.Season);

      if (!seasonMonths[seasonName]) return;

      var comp = makeComposite(year, seasonName);
      var tag = seasonName + '_' + year;

      if (EXPORT_NDVI) {
        var ndviParams = {
          image: comp.select('NDVI').toFloat(),
          description: 'Angren_NDVI_' + tag,
          folder: EXPORT_FOLDER,
          scale: MAP_SCALE,
          region: roi,
          maxPixels: 1e13,
          fileFormat: 'GeoTIFF'
        };
        if (EXPORT_CRS) ndviParams.crs = EXPORT_CRS;
        Export.image.toDrive(ndviParams);
      }

      if (EXPORT_LST) {
        var lstParams = {
          image: comp.select('LST_C').toFloat(),
          description: 'Angren_LST_' + tag,
          folder: EXPORT_FOLDER,
          scale: MAP_SCALE,
          region: roi,
          maxPixels: 1e13,
          fileFormat: 'GeoTIFF'
        };
        if (EXPORT_CRS) lstParams.crs = EXPORT_CRS;
        Export.image.toDrive(lstParams);
      }

      if (EXPORT_MASK) {
        var maskParams = {
          image: comp.select('Red').mask().rename('ValidMask').toUint8(),
          description: 'Angren_ValidMask_' + tag,
          folder: EXPORT_FOLDER,
          scale: MAP_SCALE,
          region: roi,
          maxPixels: 1e13,
          fileFormat: 'GeoTIFF'
        };
        if (EXPORT_CRS) maskParams.crs = EXPORT_CRS;
        Export.image.toDrive(maskParams);
      }
    });

    print('All export tasks created. Open the Tasks tab and run them.');
  });
}

// =============================================================================
// ANGREN: QC + PROOF EDITION (1990–2024)
// Corrected export script
// ----------------------------------------------------------------------------
// This script produces a single master table (finalData) of QC and proof
// metrics across all years and seasons, and exports quality-filtered subsets.
// It also exports per-scene QC tables and a study area pack (ROI vector,
// hillshade, DEM, static terrain zones). Finally it creates composite maps
// (NDVI, LST_C and valid pixel mask) for each year–season combination where
// enough scenes and composite coverage are available. Tasks are named
// consistently and no duplicate exports are created.
//
// IMPORTANT: Do not set the export CRS to a server-side projection. Use a
// literal CRS string (e.g., EPSG:32642) or omit the CRS completely. This
// script sets EXPORT_CRS to 'EPSG:32642' (UTM zone 42N) for all maps.
// =============================================================================

// -----------------------------------------------------------------------------
// 1) SETTINGS
// -----------------------------------------------------------------------------
var roi = ee.Geometry.Rectangle([70.1446, 41.0091, 70.2158, 41.0632]);
Map.centerObject(roi, 12);
Map.addLayer(roi, {color: 'red', fillColor: '00000000'}, 'ROI Angren');

var startYear = 1990;
var endYear   = 2024;
var years     = ee.List.sequence(startYear, endYear);

var seasons = [
  {name: 'Spring', start: 3, end: 5,  order: 1},
  {name: 'Summer', start: 6, end: 8,  order: 2},
  {name: 'Autumn', start: 9, end: 11, order: 3}
];

var MIN_PX = 10;                 // minimum valid pixels per scene
var MIN_SCENES_PROOF = 3;        // minimum scenes to compute proof metrics
var COV_THRESHOLD = 0.5;         // minimum composite coverage to compute proof metrics

// -----------------------------------------------------------------------------
// 2) REFERENCE PROJECTION (LOCK) + ROI AREA IN THAT PROJECTION
// -----------------------------------------------------------------------------
var refImg = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
  .merge(ee.ImageCollection('LANDSAT/LE07/C02/T1_L2'))
  .merge(ee.ImageCollection('LANDSAT/LC08/C02/T1_L2'))
  .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2'))
  .filterBounds(roi)
  .first();

var refProj = ee.Image(refImg).select('QA_PIXEL').projection();

var roiAreaGeo = ee.Number(roi.area(1));
var roiArea30m = ee.Number(
  ee.Image.pixelArea().reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: roi,
    scale: 30,
    crs: refProj,
    maxPixels: 1e9,
    tileScale: 4
  }).get('area')
);

// -----------------------------------------------------------------------------
// 3) HELPERS: Modern sensor check & Masking
// -----------------------------------------------------------------------------
function isL8L9(image) {
  return image.bandNames().contains('ST_B10');
}

function maskLandsatStrict(image) {
  var qa = image.select('QA_PIXEL');

  var baseMask = qa.bitwiseAnd(1 << 0).eq(0)
    .and(qa.bitwiseAnd(1 << 1).eq(0))
    .and(qa.bitwiseAnd(1 << 3).eq(0))
    .and(qa.bitwiseAnd(1 << 4).eq(0))
    .and(qa.bitwiseAnd(1 << 5).eq(0));

  var modern = isL8L9(image);
  var finalMask = ee.Image(ee.Algorithms.If(
    modern,
    baseMask.and(qa.bitwiseAnd(1 << 2).eq(0)), // cirrus only for L8/9
    baseMask
  ));

  return image.updateMask(finalMask);
}

// -----------------------------------------------------------------------------
// 4) SCALE FACTORS + BAND HARMONIZATION + INDICES
// -----------------------------------------------------------------------------
function applyScaleFactors(image) {
  var optical = image.select('SR_B.').multiply(0.0000275).add(-0.2);
  var thermal = image.select('ST_B.*').multiply(0.00341802).add(149.0);
  return image.addBands(optical, null, true)
              .addBands(thermal, null, true);
}

function unifyBands(image) {
  var modern = isL8L9(image);
  var bands = ee.List(ee.Algorithms.If(
    modern,
    ['SR_B4', 'SR_B5', 'ST_B10'],  // L8/9
    ['SR_B3', 'SR_B4', 'ST_B6']    // L5/7
  ));
  return image.select(bands, ['Red', 'NIR', 'Temp'])
    .copyProperties(image, ['system:time_start', 'LANDSAT_SCENE_ID', 'SPACECRAFT_ID']);
}

function addIndices(image) {
  var ndvi = image.normalizedDifference(['NIR', 'Red']).rename('NDVI');
  var savi = image.expression(
    '((NIR - RED) / (NIR + RED + 0.5)) * 1.5',
    {'NIR': image.select('NIR'), 'RED': image.select('Red')}
  ).rename('SAVI');
  var lstC = image.select('Temp').subtract(273.15).rename('LST_C');
  return image.addBands([ndvi, savi, lstC]);
}

// -----------------------------------------------------------------------------
// 5) TERRAIN (NASADEM) + ZONES
// -----------------------------------------------------------------------------
var dem = ee.Image('NASA/NASADEM_HGT/001');
var elevation = dem.select('elevation');
var slope = ee.Terrain.slope(elevation).rename('Slope');
var aspectRad = ee.Terrain.aspect(elevation).multiply(Math.PI).divide(180);
var northness = aspectRad.cos().rename('Northness'); // +1 north, -1 south
var terrainStack = elevation.rename('elevation')
  .addBands(slope)
  .addBands(northness);

var elevPerc = elevation.reduceRegion({
  reducer: ee.Reducer.percentile([33, 66]),
  geometry: roi,
  scale: 30,
  crs: refProj,
  maxPixels: 1e9,
  tileScale: 4
});
var p33 = ee.Number(elevPerc.get('elevation_p33'));
var p66 = ee.Number(elevPerc.get('elevation_p66'));
var zoneLow  = elevation.lt(p33);
var zoneMid  = elevation.gte(p33).and(elevation.lt(p66));
var zoneHigh = elevation.gte(p66);

// -----------------------------------------------------------------------------
// 6) ERA5-LAND (native ~11 km) aggregation per year-season
// -----------------------------------------------------------------------------
function aggregateClimate(year, season) {
  var eraCol = ee.ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR')
    .filterBounds(roi)
    .filter(ee.Filter.calendarRange(year, year, 'year'))
    .filter(ee.Filter.calendarRange(season.start, season.end, 'month'));

  var pp_m = eraCol.select('total_precipitation_sum').sum().reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: roi,
      scale: 11132,
      maxPixels: 1e9
    }).get('total_precipitation_sum');
  pp_m = ee.Algorithms.If(pp_m, pp_m, 0);
  var pp_mm = ee.Number(pp_m).multiply(1000);

  var ssm = eraCol.select('volumetric_soil_water_layer_1').mean().reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: roi,
      scale: 11132,
      maxPixels: 1e9
    }).get('volumetric_soil_water_layer_1');
  ssm = ee.Algorithms.If(ssm, ssm, 0);

  return ee.Dictionary({'PP_mm': pp_mm, 'SSM_m3m3': ssm});
}

// -----------------------------------------------------------------------------
// 7) HELPERS for robust reduceRegion means with null
// -----------------------------------------------------------------------------
function meanBand(img, bandName, maskImg) {
  var x = img.select(bandName);
  x = ee.Image(ee.Algorithms.If(maskImg, x.updateMask(maskImg), x));
  var v = x.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: roi,
    scale: 30,
    crs: refProj,
    maxPixels: 1e9,
    tileScale: 4
  }).get(bandName);
  return ee.Algorithms.If(v, v, null);
}

// -----------------------------------------------------------------------------
// 8) MAIN LOOP: YEAR x SEASON → finalData
// -----------------------------------------------------------------------------
var finalResults = ee.FeatureCollection(
  years.map(function(y) {
    return ee.List(seasons).map(function(seasonObj) {
      var season = ee.Dictionary(seasonObj);
      var year = ee.Number(y);

      var climate = aggregateClimate(year, {
        start: season.getNumber('start'),
        end: season.getNumber('end')
      });

      var col = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
        .merge(ee.ImageCollection('LANDSAT/LE07/C02/T1_L2'))
        .merge(ee.ImageCollection('LANDSAT/LC08/C02/T1_L2'))
        .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2'))
        .filterBounds(roi)
        .filter(ee.Filter.calendarRange(year, year, 'year'))
        .filter(ee.Filter.calendarRange(season.getNumber('start'), season.getNumber('end'), 'month'));

      var withPx = col.map(function(img) {
        var masked = maskLandsatStrict(img);
        var px = masked.select('QA_PIXEL').reduceRegion({
          reducer: ee.Reducer.count(),
          geometry: roi,
          scale: 30,
          crs: refProj,
          maxPixels: 1e9,
          tileScale: 4
        }).get('QA_PIXEL');
        px = ee.Number(ee.Algorithms.If(px, px, 0));
        return img.set('valid_pixels', px);
      });

      var validImages = withPx.filter(ee.Filter.gte('valid_pixels', MIN_PX));
      var count = validImages.size();

      var validWithCoverage = validImages.map(function(img) {
        var masked = maskLandsatStrict(img);
        var sceneArea = ee.Image.pixelArea().updateMask(masked.select('QA_PIXEL').mask()).reduceRegion({
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

      var covSceneMean = ee.Number(ee.Algorithms.If(
        count.gt(0),
        validWithCoverage.aggregate_mean('scene_cov_area_ratio'),
        0
      ));

      return ee.Algorithms.If(count.gt(0), (function() {
        var processed = ee.ImageCollection(validWithCoverage)
          .map(maskLandsatStrict)
          .map(applyScaleFactors)
          .map(unifyBands)
          .map(addIndices);

        var composite = processed.median().clip(roi);

        var compArea = ee.Image.pixelArea().updateMask(composite.select('Red').mask()).reduceRegion({
          reducer: ee.Reducer.sum(),
          geometry: roi,
          scale: 30,
          crs: refProj,
          maxPixels: 1e9,
          tileScale: 4
        }).get('area');
        compArea = ee.Number(ee.Algorithms.If(compArea, compArea, 0));
        var covComp = compArea.divide(roiArea30m).min(1).max(0);

        var stats = composite.select(['NDVI','SAVI','LST_C']).reduceRegion({
          reducer: ee.Reducer.mean().combine({
            reducer2: ee.Reducer.stdDev(),
            sharedInputs: true
          }),
          geometry: roi,
          scale: 30,
          crs: refProj,
          maxPixels: 1e9,
          tileScale: 4
        });

        var cond = covComp.gte(COV_THRESHOLD).and(count.gte(MIN_SCENES_PROOF));
        var analysisStack = composite.addBands(terrainStack);

        var lstGrad100m = ee.Algorithms.If(cond, (function() {
           var regression = analysisStack.select(['elevation', 'LST_C']).reduceRegion({
             reducer: ee.Reducer.linearFit(),
             geometry: roi,
             scale: 30,
             crs: refProj,
             maxPixels: 1e9,
             tileScale: 4
           });
           return ee.Number(ee.Algorithms.If(
             regression.get('scale'),
             ee.Number(regression.get('scale')).multiply(100),
             null
           ));
        })(), null);

        var maskNorth = northness.gt(0.3);
        var maskSouth = northness.lt(-0.3);

        var LST_North = ee.Algorithms.If(
          cond,
          meanBand(analysisStack, 'LST_C', maskNorth),
          null
        );
        var LST_South = ee.Algorithms.If(
          cond,
          meanBand(analysisStack, 'LST_C', maskSouth),
          null
        );
        var dLST_South_North = ee.Algorithms.If(
          cond,
          ee.Algorithms.If(
            LST_North,
            ee.Algorithms.If(LST_South, ee.Number(LST_South).subtract(ee.Number(LST_North)), null),
            null
          ),
          null
        );

        var NDVI_Low_Zone = ee.Algorithms.If(cond, meanBand(analysisStack, 'NDVI', zoneLow), null);
        var NDVI_Mid_Zone = ee.Algorithms.If(cond, meanBand(analysisStack, 'NDVI', zoneMid), null);
        var NDVI_High_Zone = ee.Algorithms.If(cond, meanBand(analysisStack, 'NDVI', zoneHigh), null);

        return ee.Feature(null, {
          'Year': year,
          'Season': season.getString('name'),
          'SeasonOrder': season.getNumber('order'),

          'Image_Count': count,
          'MIN_PX': MIN_PX,
          'ROI_Area_Geo_m2': roiAreaGeo,
          'ROI_Area_30m_m2': roiArea30m,

          'Cov_Scene_Mean': covSceneMean,
          'Cov_Composite': covComp,

          'NDVI_mean': stats.get('NDVI_mean'),
          'NDVI_std': stats.get('NDVI_stdDev'),
          'SAVI_mean': stats.get('SAVI_mean'),
          'SAVI_std': stats.get('SAVI_stdDev'),
          'LST_mean_C': stats.get('LST_C_mean'),
          'LST_std_C': stats.get('LST_C_stdDev'),

          'PP_mm': climate.get('PP_mm'),
          'SSM_m3m3': climate.get('SSM_m3m3'),

          'LST_Grad_100m': lstGrad100m,
          'LST_North': LST_North,
          'LST_South': LST_South,
          'dLST_South_North': dLST_South_North,
          'NDVI_Low_Zone': NDVI_Low_Zone,
          'NDVI_Mid_Zone': NDVI_Mid_Zone,
          'NDVI_High_Zone': NDVI_High_Zone
        });

      })(), ee.Feature(null, {
        'Year': year,
        'Season': season.getString('name'),
        'SeasonOrder': season.getNumber('order'),
        'Image_Count': 0,
        'MIN_PX': MIN_PX,
        'Cov_Scene_Mean': 0,
        'Cov_Composite': 0,
        'PP_mm': climate.get('PP_mm'),
        'SSM_m3m3': climate.get('SSM_m3m3')
      }));
    });
  }).flatten()
);

var finalData = finalResults.sort('Year').sort('SeasonOrder');

// -----------------------------------------------------------------------------
// 9) MASTER EXPORTS
// -----------------------------------------------------------------------------
// Export the full finalData table
Export.table.toDrive({
  collection: finalData,
  description: 'Angren_QC_Proof_Dataset',
  fileFormat: 'CSV',
  selectors: [
    'Year','Season','SeasonOrder','Image_Count','MIN_PX',
    'ROI_Area_Geo_m2','ROI_Area_30m_m2',
    'Cov_Scene_Mean','Cov_Composite',
    'NDVI_mean','NDVI_std',
    'SAVI_mean','SAVI_std',
    'LST_mean_C','LST_std_C',
    'PP_mm','SSM_m3m3',
    'LST_Grad_100m','LST_North','LST_South','dLST_South_North',
    'NDVI_Low_Zone','NDVI_Mid_Zone','NDVI_High_Zone'
  ]
});

// Export QC-filtered subsets (optional but recommended for analysis)
var QC80 = finalData.filter(ee.Filter.gt('Image_Count', 0)).filter(ee.Filter.gte('Cov_Composite', 0.80));
var QC90 = finalData.filter(ee.Filter.gt('Image_Count', 0)).filter(ee.Filter.gte('Cov_Composite', 0.90));
var PROOF_READY = finalData
  .filter(ee.Filter.gt('Image_Count', 0))
  .filter(ee.Filter.gte('Cov_Composite', COV_THRESHOLD))
  .filter(ee.Filter.gte('Image_Count', MIN_SCENES_PROOF));

Export.table.toDrive({ collection: QC80,        description: 'Angren_QC_Proof_Dataset_QC80',      fileFormat: 'CSV' });
Export.table.toDrive({ collection: QC90,        description: 'Angren_QC_Proof_Dataset_QC90',      fileFormat: 'CSV' });
Export.table.toDrive({ collection: PROOF_READY, description: 'Angren_QC_Proof_Dataset_ProofReady', fileFormat: 'CSV' });

// -----------------------------------------------------------------------------
// 10) SUPPLEMENTARY: PER-SCENE QC TABLE
// -----------------------------------------------------------------------------
var perSceneQC = ee.FeatureCollection(
  years.map(function(y) {
    y = ee.Number(y);
    return ee.FeatureCollection(
      ee.List(seasons).map(function(seasonObj) {
        var season = ee.Dictionary(seasonObj);
        var col = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
          .merge(ee.ImageCollection('LANDSAT/LE07/C02/T1_L2'))
          .merge(ee.ImageCollection('LANDSAT/LC08/C02/T1_L2'))
          .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2'))
          .filterBounds(roi)
          .filter(ee.Filter.calendarRange(y, y, 'year'))
          .filter(ee.Filter.calendarRange(season.getNumber('start'), season.getNumber('end'), 'month'));

        var fc = ee.FeatureCollection(col.map(function(img) {
          var masked = maskLandsatStrict(img);
          var validPx = masked.select('QA_PIXEL').reduceRegion({
            reducer: ee.Reducer.count(),
            geometry: roi,
            scale: 30,
            crs: refProj,
            maxPixels: 1e9,
            tileScale: 4
          }).get('QA_PIXEL');
          validPx = ee.Number(ee.Algorithms.If(validPx, validPx, 0));

          var sceneArea = ee.Image.pixelArea().updateMask(masked.select('QA_PIXEL').mask()).reduceRegion({
            reducer: ee.Reducer.sum(),
            geometry: roi,
            scale: 30,
            crs: refProj,
            maxPixels: 1e9,
            tileScale: 4
          }).get('area');
          sceneArea = ee.Number(ee.Algorithms.If(sceneArea, sceneArea, 0));
          var ratio = sceneArea.divide(roiArea30m).max(0).min(1);

          return ee.Feature(null, {
            'system_time_start': img.date().format('YYYY-MM-dd'),
            'Year': y,
            'Season': season.getString('name'),
            'SeasonOrder': season.getNumber('order'),
            'SPACECRAFT': img.get('SPACECRAFT_ID'),
            'SCENE_ID': img.get('LANDSAT_SCENE_ID'),
            'Valid_Pixels': validPx,
            'Coverage_Ratio': ratio
          });
        }));
        fc = fc.filter(ee.Filter.gte('Valid_Pixels', MIN_PX));
        return fc;
      })
    ).flatten();
  })
).flatten();

Export.table.toDrive({
  collection: perSceneQC,
  description: 'Angren_Supplementary_PerScene_QC',
  fileFormat: 'CSV',
  selectors: [
    'system_time_start','Year','Season','SeasonOrder',
    'SPACECRAFT','SCENE_ID','Valid_Pixels','Coverage_Ratio'
  ]
});

// -----------------------------------------------------------------------------
// 11) STUDY AREA PACK (Figure 1 context)
// -----------------------------------------------------------------------------
var EXPORT_FOLDER = 'Angren_JournalPack';
var MAP_SCALE = 30;
// Fixed CRS string for all map exports (UTM 42N). Set to null to use default.
var EXPORT_CRS = 'EPSG:32642';

if (true) { // study area vector
  var roiFc = ee.FeatureCollection([ee.Feature(roi, {'name': 'ROI_Angren'})]);
  Export.table.toDrive({
    collection: roiFc,
    description: 'Angren_StudyArea_ROI',
    folder: EXPORT_FOLDER,
    fileFormat: 'SHP'
  });
}

if (true) { // hillshade
  var hillshade = ee.Terrain.hillshade(elevation).rename('Hillshade');
  var hsParams = {
    image: hillshade.toUint8(),
    description: 'Angren_StudyArea_Hillshade',
    folder: EXPORT_FOLDER,
    scale: MAP_SCALE,
    region: roi,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF'
  };
  if (EXPORT_CRS) hsParams.crs = EXPORT_CRS;
  Export.image.toDrive(hsParams);
}

if (true) { // DEM
  var demParams = {
    image: elevation.rename('Elevation_m').toFloat(),
    description: 'Angren_StudyArea_DEM',
    folder: EXPORT_FOLDER,
    scale: MAP_SCALE,
    region: roi,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF'
  };
  if (EXPORT_CRS) demParams.crs = EXPORT_CRS;
  Export.image.toDrive(demParams);
}

if (true) { // static terrain + zones
  var staticMaps = terrainStack
    .addBands(zoneLow.rename('Zone_Low').toUint8())
    .addBands(zoneMid.rename('Zone_Mid').toUint8())
    .addBands(zoneHigh.rename('Zone_High').toUint8());
  var stParams = {
    image: staticMaps,
    description: 'Angren_Static_Terrain_Zones',
    folder: EXPORT_FOLDER,
    scale: MAP_SCALE,
    region: roi,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF'
  };
  if (EXPORT_CRS) stParams.crs = EXPORT_CRS;
  Export.image.toDrive(stParams);
}

// -----------------------------------------------------------------------------
// 12) FULL MAP EXPORTS: NDVI, LST_C, ValidMask per year × season
// -----------------------------------------------------------------------------
// Export maps only when a season-year combination exists in finalData (and
// passes QC filters). By default we require at least 3 scenes and 80% coverage.

// Map export parameters
var MAPS_MIN_SCENES = 3;
var MAPS_MIN_COV    = 0.80;

function makeComposite(year, seasonName) {
  var m = (seasonName === 'Spring') ? {start: 3, end: 5} :
          (seasonName === 'Summer') ? {start: 6, end: 8} : {start: 9, end: 11};
  var col = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
    .merge(ee.ImageCollection('LANDSAT/LE07/C02/T1_L2'))
    .merge(ee.ImageCollection('LANDSAT/LC08/C02/T1_L2'))
    .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2'))
    .filterBounds(roi)
    .filter(ee.Filter.calendarRange(year, year, 'year'))
    .filter(ee.Filter.calendarRange(m.start, m.end, 'month'))
    .map(maskLandsatStrict)
    .map(applyScaleFactors)
    .map(unifyBands)
    .map(addIndices);
  return col.median().clip(roi);
}

function exportAllYearSeasonMapsFromFinalData() {
  var fc = finalData
    .filter(ee.Filter.gt('Image_Count', 0))
    .filter(ee.Filter.gte('Image_Count', MAPS_MIN_SCENES))
    .filter(ee.Filter.gte('Cov_Composite', MAPS_MIN_COV))
    .select(['Year', 'Season']);

  fc.toList(fc.size()).evaluate(function(list) {
    if (!list || list.length === 0) {
      print('No year-season combinations passed the export filters.');
      return;
    }
    print('Creating map export tasks for', list.length, 'year-season combinations');
    list.forEach(function(f) {
      var p = f.properties;
      var year = Number(p.Year);
      var seasonName = String(p.Season);
      var comp = makeComposite(year, seasonName);
      var tag = seasonName + '_' + year;
      // NDVI
      var ndviParams = {
        image: comp.select('NDVI').toFloat(),
        description: 'Angren_NDVI_' + tag,
        folder: EXPORT_FOLDER,
        scale: MAP_SCALE,
        region: roi,
        maxPixels: 1e13,
        fileFormat: 'GeoTIFF'
      };
      if (EXPORT_CRS) ndviParams.crs = EXPORT_CRS;
      Export.image.toDrive(ndviParams);
      // LST
      var lstParams = {
        image: comp.select('LST_C').toFloat(),
        description: 'Angren_LST_' + tag,
        folder: EXPORT_FOLDER,
        scale: MAP_SCALE,
        region: roi,
        maxPixels: 1e13,
        fileFormat: 'GeoTIFF'
      };
      if (EXPORT_CRS) lstParams.crs = EXPORT_CRS;
      Export.image.toDrive(lstParams);
      // ValidMask
      var maskParams = {
        image: comp.select('Red').mask().rename('ValidMask').toUint8(),
        description: 'Angren_ValidMask_' + tag,
        folder: EXPORT_FOLDER,
        scale: MAP_SCALE,
        region: roi,
        maxPixels: 1e13,
        fileFormat: 'GeoTIFF'
      };
      if (EXPORT_CRS) maskParams.crs = EXPORT_CRS;
      Export.image.toDrive(maskParams);
    });
    print('All map export tasks created. Run them from the Tasks tab.');
  });
}

// Uncomment the next line to enqueue map exports for all valid year-season combos
exportAllYearSeasonMapsFromFinalData();

