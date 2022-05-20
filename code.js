
// Loading Area of Interest ---- AOI-----Area of study
var AOI = geometry;
var newfc = table;
 Map.addLayer (newfc);
//print('ref', RefData);

// Load the Sentinel 2 MSI reflectance Image Collection.
var sentinelCollection = imageCollection.filterDate('2018-01-01', '2018-12-31');
################################################################################

// Make a cloud-free composite.
// Function to mask clouds using the Sentinel-2 QA band.
function maskS2clouds(image) {
  var qa = image.select('QA60');
  // Bits 10 and 11 are clouds and cirrus, respectively.
  var cloudBitMask = ee.Number(2).pow(10).int();
  var cirrusBitMask = ee.Number(2).pow(11).int();
  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0).and(
             qa.bitwiseAnd(cirrusBitMask).eq(0));
  // Return the masked and scaled data.
  return image.updateMask(mask).divide(10000);}
  
  // Map the function during the long rain season of data and take the median.
var composite = sentinelCollection
                  // Pre-filter to get less cloudy granules.
                  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
                  .filterBounds(AOI)
                  .map(maskS2clouds)
                  .median();
                  
//print('composite', composite)

// Display the composite
//Map.addLayer(composite.clip(AOI), {bands: ['B11', 'B8', 'B4'], min: 0, max: 0.3}, 'composite');
################################################################################

//VEGETATION INDICES COMPUTATION
      //NDVI
  var NDVI = composite.normalizedDifference(['B8', 'B4'])
  
      //EVI
 var EVI = composite.expression(
    '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
    'NIR':composite.select('B8'),
    'RED':composite.select('B4'),
    'BLUE':composite.select('B2')
  }); 
  
      //MSAVI
     var MSAVI = composite.expression(
  '(2 * NIR + 1 - sqrt(pow((2 * NIR + 1), 2) - 8 * (NIR - RED)) ) / 2', {
    'NIR':composite.select('B8'), 
    'RED':composite.select('B4')

//RED-EDGE VEGETATION INDICES COMPUTATION
//Resampling 
//(https://github.com/sacridini/GEET)
// (https://github.com/sacridini/GEET/blob/master/geet.js)
var resample_band = function (band, scale, mode) {
    // Error Handling
    //if (image === undefined) error('resample_band', 'You need to specify an input image.');
    //if (scale === undefined) error('resample_band', 'You need to specify the scale number.');
	//if (mode === undefined) error('resample', 'You need to specify the resample mode (bilinear or bicubic).');
    var resampled_B6 = band.resample(bilinear).reproject({
        crs: band.projection().crs(),
        scale:10
    });
    return resampled_band;
};

// Write here functions that add bands to the image collection or carry out some kind of processing e.g. as shown below
// Function to calculate and add msavi band. This will be executed in the map_and_export_EVI function below. 

//Red_Edge_EVI

var EVI_RE = composite.expression(
    '2.5 * ((NIR - RED_EDGE) / (NIR + 6 * RED_EDGE - 7.5 * BLUE + 1))', {
    'NIR':composite.select('B8'),
    'RED_EDGE':composite.select('B6'),
    'BLUE':composite.select('B2')
  }); 

//Red_Edge_NDVI 
var NDVI_RE = composite.normalizedDifference(['B8', 'B6']);
################################################################################

//VEGETATION PHENOMETRICS DEC-FEB 2018-2019
var phenology = image;
//print('phenology', phenology);
var b1 = phenology.select('b1')       //Onset_Value (Gray)
var b2 = phenology.select('b2')       //Onset_Time
var b3 = phenology.select('b3')       //Max_Value
var b4 = phenology.select('b4')       //Max_Time
var b5 = phenology.select('b5')       //Offset_Value
var b6 = phenology.select('b6')       //Offse_Time
var b7 = phenology.select('b7')       //Length_GS
var b8 = phenology.select('b8')       //Before Max T
var b9 = phenology.select('b9')       //After Max T
var b10 = phenology.select('b10')     //Green Up Space
var b11 = phenology.select('b11')     //Brown Down Slope
var b12 = phenology.select('b12')     //TINDVI Before Max
var b13 = phenology.select('b13')     //TINDVI After Max        
var b14 = phenology.select('b14')     //TINDVI
var b15 = phenology.select('b15')     //Asymmetry
################################################################################

//COMBINATION OF ALL THE BANDS AND IMAGES
var composite2 = composite.addBands(NDVI.clip(AOI).rename('ndvi'))
                   .addBands(EVI.clip(AOI).rename('evi'))
                   .addBands(MSAVI.clip(AOI).rename('msavi'))
                   .addBands(EVI_RE.clip(AOI).rename('evi_re'))
                   .addBands(NDVI_RE.clip(AOI).rename('ndvi_re'))
                    .addBands (b1.clip(AOI).rename('b1_pheno'))
                    .addBands (b2.clip(AOI).rename('b2_pheno'))
                    .addBands (b3.clip(AOI).rename('b3_pheno'))
                    .addBands (b4.clip(AOI).rename('b4_pheno'))
                    .addBands (b5.clip(AOI).rename('b5_pheno'))
                    .addBands (b6.clip(AOI).rename('b6_pheno'))
                    .addBands (b7.clip(AOI).rename('b7_pheno'))
                    .addBands (b8.clip(AOI).rename('b8_pheno'))
                    .addBands (b9.clip(AOI).rename('b9_pheno'))
                    .addBands (b10.clip(AOI).rename('b10_pheno'))
                    .addBands (b11.clip(AOI).rename('b11_pheno'))
                    .addBands (b12.clip(AOI).rename('b12_pheno'))
                    .addBands (b13.clip(AOI).rename('b13_pheno'))
                    .addBands (b14.clip(AOI).rename('b14_pheno'))
                    .addBands (b15.clip(AOI).rename('b15_pheno'))
   //print('composite2', composite2);
################################################################################
Mapping Avocado Agro-ecological Systems using Sentinel-2
//LULC CLASSIFICATION 
// Merge the five geometry layers into a single FeatureCollection.
//Fusion table in use. Point data. learn more from here (https://mygeoblog.com/2016/10/17/creating-fusion-tables/)
//Fusion table in use. Polygon data. learn more from here (https://https://www.gislounge.com/converting-shapefiles-to-google-fusion-tables/)

//var newfc = RefData;
//var newfc = RefData.merge(Water);

// Use these bands for classification.
//VI inclusive
//var bands = ['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A','evi','msavi','ndvi'];

//RED_EDGE Inclusive
//var bands = ['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A','evi_re','ndvi_re'];

//VI & RED_EDGE
//var bands = ['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A','evi','msavi','ndvi','evi_re','ndvi_re'];

/*
//Pheometrics inclusive
var bands = ['B2',
              'B3', 
              'B4', 
              'B5', 
              'B6', 
              'B7', 
              'B8', 
              'B8A',
              'b1_pheno',
              'b2_pheno',
              'b3_pheno',
              'b4_pheno',
              'b5_pheno',
              'b6_pheno',
              'b7_pheno',
              'b8_pheno',
              'b9_pheno',
              'b10_pheno',
              'b11_pheno',
              'b12_pheno',
              'b13_pheno',
              'b14_pheno',
              'b15_pheno',];
*/

//VI + Red_Edge + Phenology inclusive
var bands = ['B2',
              'B3', 
              'B4', 
              'B5', 
              'B6', 
              'B7', 
              'B8', 
              'B8A',
              'evi',
              'msavi',
              'ndvi',
              'evi_re',
              'ndvi_re',
              'b1_pheno',
              'b2_pheno',
              'b3_pheno',
              'b4_pheno',
              'b5_pheno',
              'b6_pheno',
              'b7_pheno',
              'b8_pheno',
              'b9_pheno',
              'b10_pheno',
              'b11_pheno',
              'b12_pheno',
              'b13_pheno',
              'b14_pheno',
              'b15_pheno',];

/*
// Red_Edge + Phenology Band (3,4,7,12,13,14,15) inclusive
var bands = ['B2',
              'B3', 
              'B4', 
              'B5', 
              'B6', 
              'B7', 
              'B8', 
              'B8A',
              'evi_re',
              'ndvi_re',
              'b3_pheno',
              'b4_pheno',
              'b7_pheno',
              'b12_pheno',
              'b13_pheno',
              'b14_pheno',
              'b15_pheno',];
*/

// Split into training and testing features
// The name of the property on the points storing the class label.
var classProperty = 'Landcover';

// Optionally, do some accuracy assessment.  Fist, add a column of
// random uniforms to the training dataset.
var withRandom = newfc.randomColumn('random');
// We want to reserve some of the data for testing, to avoid overfitting the model.
var split = 0.7;  // Roughly 70% training, 30% testing.
//var trainingPartition = withRandom.filter(ee.Filter.lt('random', split));
//var testingPartition = withRandom.filter(ee.Filter.gte('random', split));

var trainingPartition = withRandom.filter(ee.Filter.lt('random', split));
var testingPartition = withRandom.filter(ee.Filter.gte('random', split));
//print('trinP',trainingPartition)
//print('testP',testingPartition)

// Sample the composite to generate training data.  Note that the
// class label is stored in the 'landcover' property.
 
 var training = composite2.select(bands).sampleRegions({
  collection: trainingPartition,
  properties: [classProperty],
  scale: 10
});
print('train', training);

// Train a CART classifier.//chnage ro Random Forest 
var classifier = ee.Classifier.smileRandomForest({
            numberOfTrees: 100,//was less
            seed: 1
          })
        .train({
            features: training,
            classProperty: classProperty,
});

// Classify the composite.
var classified = composite2.classify(classifier);
var palette = [
  'ff78cb', // Annua_Croplands     (0) Pink
  '6b8624', // Avocado             (1) olive
  'ff0000', // Built_Up            (2) red
  'FFFF00', // Grasslands          (3) yellow
  'ff7f50', // Perrenial_Croplands (4) coral
  '8b4513', // Shrublands          (5) brown
  '008b00', // Tree_cover          (6) green
  '0000FF' //  water               (7) blue
];
Map.addLayer(classified.clip(AOI), {min: 1, max: 8, palette: palette}, 'Land Use Classification');
Map.centerObject(newfc);

 
 var testing = composite2.select(bands).sampleRegions({
  collection: testingPartition,
  properties: [classProperty],
  scale: 10
});
print('testing', testing)

// Classify the test FeatureCollection.
var test = testing.classify(classifier);


// Compute pixel area in square meters per landcover type.
var stats = ee.Image.pixelArea().addBands(classified).reduceRegion({
          reducer: ee.Reducer.sum().group(1),
          geometry:AOI,
          maxPixels: 1e13,
          tileScale:2,
          scale: 10,
          //crs: "EPSG:4326"
        });
 
################################################################################
//VALIDATION & ACCURACY
var confMatrix = classifier.confusionMatrix();
//print(confMatrix);

var testAccuracy = test.errorMatrix('Landcover', 'classification');
var OA = confMatrix.accuracy();
var CA = confMatrix.consumersAccuracy();
var Kappa = confMatrix.kappa();
var Order = confMatrix.order();
var PA = confMatrix.producersAccuracy();

//print(confMatrix,'Confusion Matrix');
//print(OA,'Overall Accuracy');
//print(CA,'Consumers Accuracy');
//print(Kappa,'Kappa');
//print(Order,'Order');
//print(PA,'Producers Accuracy');

print('Test error matrix: ', testAccuracy);
print('Test overall accuracy: ', testAccuracy.accuracy());
print('Test Kappa: ', testAccuracy.kappa());
print('Test consumersAccuracy: ', testAccuracy.consumersAccuracy());
print('Test producersAccuracy: ', testAccuracy.producersAccuracy());
 
//print('Area per class', stats);
################################################################################

//Exporting classifieed images - When satisfied with the accuracy.

//Export.image(). Though this function take several optional parameters, it is valuable to familiarize yourself with these: 
//maxPixels \96 This restricts the numbers of pixels in the exported image. By default, this value is set to 100,000,000 pixels. 
//            You can set this argument to raise or lower the limit. \931e13\94 is 10 to the 13th power (1013) and the most GEE will let you export. 
//region \96 By default, the viewport of the Code Editor is exported but you can also specify a geometry to change the export extent. 
//crs \96 The coordinate reference system for the output image. This is specified using the EPSG code. You can look up the EPSG code for your desired spatial projections at http://spatialreference.org. 
//scale \96 The resolution in meters per pixel. The native resolution for the canopy height data set is 30 arc-seconds or approximately one kilometer.
//Export the image
//the scale change according to the image resolution

/*
//REFLECTANCE
//exporting classifieed images
Export.image.toDrive({
image: classified,
description: 'Reflectance_Point_Muranga_2018',
scale: 10,
region: geometry,
maxPixels: 3E10
});
*/

/*
//REFLECTANCE + VI
//exporting classifieed images
Export.image.toDrive({
image: classified,
description: 'Reflectance_VI_Point_Muranga_2018',
scale: 10,
region: geometry,
maxPixels: 3E10
});
*/

/*
//REFLECTANCE + RED_EDGE
//exporting classifieed images
Export.image.toDrive({
image: classified,
description: 'Reflectance_Red_Edge_Point_Muranga_2018',
scale: 10,
region: geometry,
maxPixels: 3E10
});
*/

/*
//REFLECTANCE + PHENOLOGY
//exporting classifieed images
Export.image.toDrive({
image: classified,
description: 'Reflectance_Phenology_Point_Muranga_2018',
scale: 10,
region: geometry,
maxPixels: 3E10
});
*/

/*
//REFLECTANCE + VI + RED-EDGE + PHENOLOGY
//exporting classifieed images
Export.image.toDrive({
image: classified,
description: 'Reflectance_VI_RE_Phenology_Point_Muranga_2018',
scale: 10,
region: geometry,
maxPixels: 3E10
});
*/

/*
//REFLECTANCE + VI + RED-EDGE + PHENOLOGY bands( 3, 4, 7, 12,13,14,15)
//exporting classifieed images
Export.image.toDrive({
image: classified,
description: 'Reflectance_RE_Pheno(selected)_Point_Muranga_2018',
scale: 10,
region: geometry,
maxPixels: 3E10
});
*/

/*
//REFLECTANCE + VI + RED-EDGE + PHENOLOGY bands( 3, 4, 7, 12,13,14,15)
//exporting classifieed images
Export.image.toDrive({
image: classified,
description: 'Reflectance_VI_Rededge_Point_Muranga_2018',
scale: 10,
region: geometry,
maxPixels: 3E10
});
*/
################################################################################
// Variable importance
var dict = classifier.explain();
print('Explain:',dict);

var variable_importance = ee.Feature(null, ee.Dictionary(dict).get('importance'));

var chart =
  ui.Chart.feature.byProperty(variable_importance)
    .setChartType('ColumnChart')
    .setOptions({
      title: 'All poly Random Forest Variable Importance',
      legend: {position: 'none'},
      hAxis: {title: 'Variables'},
      vAxis: {title: 'Importance'}
    });
print(chart); 
################################################################################
//Exporting Feature classes

// Export the FeatureCollection.
Export.table.toDrive({
  collection: testingPartition,
  description: 'testP',
  fileFormat: 'CSV'
});

// Export the FeatureCollection to a KML file.
Export.table.toDrive({
  collection: testingPartition,
  description:'testP',
  fileFormat: 'KML'
});

// Export the FeatureCollection to a KML file.
Export.table.toDrive({
  collection: testingPartition,
  description:'testP',
  fileFormat: 'SHP'
});


