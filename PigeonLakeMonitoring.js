//---------------------------------------------------------------------------------------------------
//Version 1 (beta): Pigeon lake chlorophyll-A prediction with Sentinel-2 and -3
///App development: Evan R. DeLancey
///email: edelance@ualberta.ca
///Spatial Data Scientist - ABMI
///
///Other team members: Caleb Sinn, Bradley Peter, Rolf Vinebrooke, Arin MacFarlane-Dyer
///Date: 2021-03-10
//-----------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------
//Map settings
//-----------------------------------------------------------------------------------------------------
Map.setControlVisibility(false, true, false, false, false, false);
var d = Map.drawingTools();
d.setShown(false);
//-----------------------------------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------------------
//Define functions
//-----------------------------------------------------------------------------------------------------
function addIndices(i){
  var a = i.addBands(i.normalizedDifference(['B5', 'B4']).rename('NDCI'));
  //threee band Ansper et al 2018
  var b = a.addBands(a.expression(
    '(B5/10000) - (((B4/10000) + (B6/10000))/2)',{
      'B4' : i.select(['B4']),
      'B5' : i.select(['B5']),
      'B6' : i.select(['B6'])
    }).rename('AnsperTBI'));
  var c = b.addBands(b.expression(
    'B5/B4',{
      'B4' : i.select(['B4']),
      'B5' : i.select(['B5']),
    }).rename('AnsperB5R'));
  var d = c.addBands(c.expression(
    'B6/B4',{
      'B4' : i.select(['B4']),
      'B6' : i.select(['B6']),
    }).rename('AnsperB6R'));
  var e = d.addBands(d.expression(
    '(AnsperTBI * 2445.115) + 5.046',{
      'AnsperTBI' : d.select(['AnsperTBI'])
    }).rename('chlA').clamp(0,300));
  return e;
}
//make image collection from slected dates S2
function makeImgCollection(f) {
  var date1 = ee.Date(f);
  var date2 = date1.advance(1, 'day');
  var img = S2.filterBounds(SA)
              .filterDate(date1, date2)
              .median();
  var img = img.clip(SA);
  var img=  img.set('date', f)
  return img.set("system:time_start", ee.Date(f).millis());
}

//make image collection from slected dates S3
function makeImgCollectionS3(f) {
  var date1 = ee.Date(f);
  var date2 = date1.advance(1, 'day');
  var img = S3.filterBounds(SA)
              .filterDate(date1, date2)
              .median();
  var img = img.clip(SAinner);
  var img=  img.set('date', f)
  return img.set("system:time_start", ee.Date(f).millis());
}

//cholorphyll A from S3
function calChlA(i){
  var MCI = i.expression(
  'L709 - L681 - 0.389*(L753 - L681)', {
    'L709': i.select(['Oa11_radiance']),
    'L681': i.select(['Oa10_radiance']),
    'L753': i.select(['Oa12_radiance'])
  }).rename(['MCI']);
  var MCIadd = i.addBands(MCI);
  var chlA = MCIadd.select('MCI').multiply(5.3109).add(6.7946);
  return i.addBands(chlA.rename('chlA').clamp(0,300));
}
//----------------------------------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------------------
//Study area and packages
//-----------------------------------------------------------------------------------------------------
var animation = require('users/gena/packages:animation');
var SA = ee.FeatureCollection('users/edelance/PigeonLake');
var colorbrewer = require('users/gena/packages:colorbrewer');
var SAinner = ee.FeatureCollection('users/edelance/PigeonLake_inner300');
//----------------------------------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------
//Map layout
//-----------------------------------------------------------------------------------------------------
Map.setOptions('HYBRID');
Map.centerObject(SA, 12);
//----------------------------------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------------------
//Pull Sentinel-2 data and filter only sunny days
//-----------------------------------------------------------------------------------------------------
var S2 = ee.ImageCollection('COPERNICUS/S2')
    .filterDate('2017-05-01', '2020-12-31')
    .filterBounds(SA)
    .map(addIndices);
//List of cloud free S2 dates 
var dates = ee.List(['2017-06-06', '2017-06-19', '2017-08-25', '2017-09-02', '2017-09-07', '2017-09-27', 
                    '2017-09-29', '2017-10-04', '2017-10-24',
                    '2018-05-15', '2018-05-20', '2018-05-22', '2018-05-25', '2018-07-16',
                    '2018-07-29', '2018-09-07', '2018-09-19', '2018-10-17', '2018-10-19', 
                    '2018-10-22', '2018-10-24', '2019-05-12', '2019-06-16', '2019-06-26', '2019-07-21',
                    '2019-07-26', '2019-08-08', '2019-08-15', '2019-08-20', '2019-08-28', '2019-09-02',
                    '2019-09-04', '2019-09-12', '2019-09-19', '2019-09-22', '2019-10-02', '2019-10-19',
                    '2020-05-24', '2020-06-20', '2020-07-10', '2020-07-28', '2020-07-30', '2020-08-02', 
                    '2020-08-04', '2020-08-09', '2020-08-19', '2020-08-27', '2020-09-11', '2020-10-01', 
                    '2020-10-03']); 
//Make image stack                    
var S2stack = ee.ImageCollection(dates.map(makeImgCollection));
//----------------------------------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------------------
//Pull Sentinel-3 and filter dates
//-----------------------------------------------------------------------------------------------------
var S3 = ee.ImageCollection('COPERNICUS/S3/OLCI')
  .filterBounds(SA)
  .filterDate('2017-05-01', '2020-10-31')
  .filterMetadata('spacecraft', 'equals', 'S3A')
  .map(calChlA);
//list of S3 dates
var dates = ee.List(['2017-05-28','2017-06-24','2017-07-02',
                    '2017-07-13','2017-07-28','2017-08-25','2017-08-28','2017-08-29','2017-09-02',
                    '2017-09-05','2017-09-06','2017-09-10','2017-09-16','2017-09-24','2017-09-28','2017-09-29',
                    '2017-10-03','2017-10-21','2017-10-22','2018-05-14','2018-05-15','2018-05-19',
                    '2018-05-22','2018-05-25','2018-05-26','2018-06-03','2018-06-18','2018-07-15',
                    '2018-07-18','2018-07-27','2018-09-06','2018-09-19','2018-09-25','2018-10-03','2018-10-11',
                    '2018-10-16','2019-05-20',
                    '2019-07-13','2019-07-23','2019-08-04','2019-09-05','2019-09-07','2019-09-12','2019-09-19',
                    '2019-10-01','2018-07-06','2018-07-07','2018-07-11','2019-10-20',
                    '2018-07-19','2018-07-26','2019-10-21','2018-07-29','2018-08-02','2018-08-06',
                    '2019-10-25','2018-09-07','2019-10-29',
                    '2020-05-13','2020-05-24','2020-06-20','2020-06-24','2020-08-01',
                    '2020-08-04','2020-08-05','2020-08-09','2019-05-27','2020-08-16','2019-07-17',
                    '2020-08-17','2020-08-23','2020-08-27','2020-08-28','2020-09-09','2020-09-16',
                    '2020-09-20','2019-10-06','2020-09-25','2020-09-27','2020-10-01',
                    '2020-10-02']); 
                 
var S3stack = ee.ImageCollection(dates.map(makeImgCollectionS3));
var S3stack = S3stack.sort('system:time_start', true);
//-----------------------------------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------------------
//Front end panels for selecting dates
//-----------------------------------------------------------------------------------------------------
var C1 = '#ffffff90';//white transparency
var C2 = '#00000000';//black full transparency
var C3 = '#000000';//black
var C7 = 'white'


var panelTitle = ui.Label({
  value: 'Selection Panel',
  style: {
    backgroundColor: C2,
    fontWeight: 'bold',
    fontSize: '18px',
  }
});

var textStartDate = ui.Textbox({
  value: '2020-05-01',
  onChange: function(){},
  style: {
    backgroundColor: 'black',
  }
});
var titleStartDate = ui.Label({
  value: 'Start date (Year-Month-date)',
  style: {
    backgroundColor: C2,
    fontWeight: 'bold',
    fontSize: '12px',
  }
});
var textEndDate = ui.Textbox({
  value: '2020-10-31',
  onChange: function(){},
  style: {
    backgroundColor: 'black',
  }
});
var titleEndDate = ui.Label({
  value: 'End date',
  style: {
    backgroundColor: C2,
    fontWeight: 'bold',
    fontSize: '12px'
  }
});
var selectBox = ui.Select({
  items: ['Chlorophyll-A', 'True colour'],
  placeholder: 'Map output',
  onChange: function() {
  }
});

var selectBoxSensor = ui.Select({
  items: ['Sentinel-2', 'Sentinel-3'],
  placeholder: 'Satellite',
  onChange: function() {
  }
});
//-----------------------------------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------------------
//Viz parameters
//-----------------------------------------------------------------------------------------------------
var pal = colorbrewer.Palettes.Spectral[11];
var pal = pal.reverse();
var vis = {bands:['B4','B3','B2'], gamma: 1.4, min:0.05, max:[2500, 2500, 3000]};
var visS3 = {bands: ['Oa08_radiance', 'Oa06_radiance', 'Oa04_radiance'], min:0, max: 100};
var visChlA = {bands:['chlA'], min:0, max: 60, palette: pal};
var visLegend = {min:0, max: 60, palette: pal};

//-----------------------------------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------------------------------




//-----------------------------------------------------------------------------------------------------
//create legend
//-----------------------------------------------------------------------------------------------------
var long = ee.Image.pixelLonLat().select('longitude');
var gradient = long.multiply((visLegend.max-visLegend.min)/100.0).add(visLegend.min);
var legendImage = gradient.visualize(visLegend);

var thumb = ui.Thumbnail({
    image: legendImage, 
    params: {bbox:'0,0,100,8', dimensions:'160x20'}, 
    style: {padding: '0px', position: 'top-center'}
});
var LC1 = '#00000050';//black transparecy
var LC2 = '#00000001';//black transparecy
var titleStyle = {
    stretch: 'vertical',
    color: 'white',
    backgroundColor: LC2,
    fontSize: '20px',
    fontWeight: 'bold'
  }
var title = ui.Panel({
    widgets: [
      ui.Label({value:'Chlorophyll A (µg/L)', style: titleStyle})
    ],
    layout: ui.Panel.Layout.flow('horizontal'),
    style: {stretch: 'horizontal', backgroundColor: LC2},
  });

var panelStyle = {
    color: 'white',
    backgroundColor: LC2,
    fontSize: '16px'
};
var panel = ui.Panel({
    widgets: [
      ui.Label({value: '0', style: panelStyle}), 
      ui.Label({style: {stretch: 'horizontal'}}),
      ui.Label({value: '30', style: panelStyle}), 
      ui.Label({style: {stretch: 'horizontal'}}),
      ui.Label({value: '> 60', style: panelStyle})
    ],
    layout: ui.Panel.Layout.flow('horizontal'),
    style: {stretch: 'horizontal', backgroundColor: LC2}
});

var legendStyle = {
    backgroundColor: C2,
    position: 'bottom-left',
    color: LC2
};
var legend = ui.Panel({style: legendStyle}).add(title).add(panel).add(thumb);
//-----------------------------------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------
//About and intructions panel from Justin Braaten EO time series explorer
//-----------------------------------------------------------------------------------------------------
var CONTROL_PANEL_WIDTH = '280px';
var CONTROL_PANEL_WIDTH_HIDE = '165px';
var infoFont = {fontSize: '11px', color: C3, backgroundColor: C2};
var titleFont = {fontSize: '13px', fontWeight: 'bold', color: C3, backgroundColor: C2};

var instrucButton = ui.Button(
  {label: 'Instructions ❯', style: {margin: '0px 4px 0px 0px'}});

var instrucTitle = ui.Label(
  'Instructions',
  titleFont);
var instrucLabel1 = ui.Label(
  '1. In the Selection Panel, enter date range (2017-present day)',
  infoFont);
var instrucLabel2 = ui.Label(
  '2. Map output will control whether users want to view chlorophyll A predictions or true colour RGB images',
  infoFont);
var instrucLabel3 = ui.Label(
  '3. The Satellite tab allows users to select Sentinel-2 or Sentinel-3 data. ' + 
  'Sentinel-2 data is 20 m resolution but less frequent images while Sentinel-3 is 300 m resolution but more frequent.',
  infoFont);
var instrucLabel4 = ui.Label(
  '4. Press “GO” to generate the images',
  infoFont);
var instrucLabel5 = ui.Label(
  '5. Once the images have loaded (the progress will be shown in the top right layers tab), you can ' + 
  'click the animation bar and use arrow keys to scroll through the dates',
  infoFont);
var instrucLabel6 = ui.Label(
  '6. A time-series chart will pop up for the desired date range.' + 
  'To download this data click the box at the top right of the chart and then select export to csv or png.',
  infoFont);

var aboutTitle = ui.Label(
  'About',
  titleFont);
var aboutButton = ui.Button(
  {label: 'About ❯', style: {margin: '0px 4px 0px 0px'}});
var aboutLabel = ui.Label(
  'This web tool helps users visualize and download Sentinel-2 and Sentinel-3 derived chlorophyll A ' + 
  'predictions for Pigeon Lake, Alberta. It is a joint project between: the Alberta Biodiversity Monitoring ' +
  'Institute, University of Alberta, Alberta Lake Management Society, Alberta Environment and Parks, and ' +
  'Pigeon Lake Watershed Association. Chlorophyll A predictions are calibrated with fields samples ' +
  'obtained in spring/summer 2020. Contact edelance@ualberta.ca for questions.',
  infoFont);
  
  
var controlPanel = ui.Panel({
  style: {position: 'bottom-right', width: CONTROL_PANEL_WIDTH_HIDE,
    maxHeight: '90%'
  }});

var buttonPanel = ui.Panel(
  [instrucButton, aboutButton],
  ui.Panel.Layout.Flow('horizontal'),
  {stretch: 'horizontal', margin: '0px 0px 0px 0px'});
  
var infoElements = ui.Panel(
  {style: {shown: false, margin: '0px -8px 0px -8px'}});
  
var aboutElements = ui.Panel(
  {style: {shown: false, margin: '0px -8px 0px -8px'}});
  
var infoShow = false;
var controlShow = false;
function infoButtonHandler() {
  if(infoShow) {
    infoShow = false;
    infoElements.style().set('shown', false);
    instrucButton.setLabel('Instructions ❯');
  } else {
    infoShow = true;
    infoElements.style().set('shown', true);
    instrucButton.setLabel('Instructions ❮');
  }
  if(infoShow | controlShow) {
    controlPanel.style().set('width', CONTROL_PANEL_WIDTH);
  } else {
    controlPanel.style().set('width', CONTROL_PANEL_WIDTH_HIDE);
  }
}

function aboutButtonHandler() {
  if(controlShow) {
    controlShow = false;
    aboutElements.style().set('shown', false);
    aboutButton.setLabel('About ❯');
  } else {
    controlShow = true;
    aboutElements.style().set('shown', true);
    aboutButton.setLabel('About ❮');
  }
  if(infoShow | controlShow) {
    controlPanel.style().set('width', CONTROL_PANEL_WIDTH);
  } else {
    controlPanel.style().set('width', CONTROL_PANEL_WIDTH_HIDE);
  }
}

infoElements.add(instrucTitle);
infoElements.add(instrucLabel1);
infoElements.add(instrucLabel2);
infoElements.add(instrucLabel3);
infoElements.add(instrucLabel4);
infoElements.add(instrucLabel5);
infoElements.add(instrucLabel6);

aboutElements.add(aboutTitle);
aboutElements.add(aboutLabel);



controlPanel.add(buttonPanel);
controlPanel.add(infoElements);
controlPanel.add(aboutElements);


Map.add(controlPanel);
instrucButton.onClick(infoButtonHandler);
aboutButton.onClick(aboutButtonHandler);
//-----------------------------------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------------------
//chl A chart options
//-----------------------------------------------------------------------------------------------------
var cw = 450;
var ch = 300;
var lw = 3;
var chartOptions = {
      interpolateNulls: true,
      vAxis: {
        title: ' Mean chlorophyll A (µg/L)', 
        titleTextStyle: {color: C3, italic: false, bold: true, fontSize:14}, 
        gridlines: {color: '#E6E6E6'}, 
        textStyle: {color: C3, italic: false, fontSize:12.5}
      },
      hAxis: {
        title: 'Day of year', 
        titleTextStyle: {color: C3, italic: false, bold: true, fontSize: 14}, 
        gridlines: {color: '#E6E6E6'}, 
        textStyle: {color: C3, italic: false, fontSize: 12.5},
    },
      series: {
        0: {color: '#B42F32'},
        1: {color: '#2D8E87'},
        2: {color: '#E8A631'},
        3: {color: '#50C878'}
      },
      backgroundColor: {fill: C7},
      curveType: 'function',
      height: ch,
      width: cw,
      fontName: 'Roboto',
      legend: {
        textStyle: {color: C3, italic: false, fontSize: 11.5}
      },
      pointSize: 4,
      lineWidth: lw,
};
//-----------------------------------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------------------
//Selector Panel function
//-----------------------------------------------------------------------------------------------------
var button = ui.Button({
  label: 'GO',
  onClick: function () {
    Map.clear();
    Map.setControlVisibility(false, true, false, false, false, false);
    var d = Map.drawingTools();
    d.setShown(false);
    //Animation package from Gennadii Donchyts
    var animation = require('users/gena/packages:animation');
    var start = String(textStartDate.getValue());
    var end = String(textEndDate.getValue());
    var output = selectBox.getValue();
    var satellite = selectBoxSensor.getValue();
    var S2filter = S2stack.filterDate(start, end);
    var S3filter = S3stack.filterDate(start, end);
    
    //get all image collections
    var imagesS2 = S2filter.map(function(i) {
      return i.visualize(vis).set({ label: i.get('date')});
    });
    var imagesS2chlA = S2filter.map(function(i) {
      return i.visualize(visChlA).set({ label: i.get('date')});
    });
    var imagesS3 = S3filter.map(function(i) {
      return i.visualize(visS3).set({ label: i.date().format('YYYY-MM-dd') });
    });
    var imagesS3chlA = S3filter.map(function(i) {
      return i.visualize(visChlA).set({ label: i.date().format('YYYY-MM-dd') });
    });
    
    //select image collection to use based on select boxes
    var images = ee.Algorithms.If(
      ee.Algorithms.IsEqual(satellite, 'Sentinel-2'),
      ee.Algorithms.If(
        ee.Algorithms.IsEqual(output, 'Chlorophyll-A'),
        imagesS2chlA,
        imagesS2
        ),
     ee.Algorithms.If(
        ee.Algorithms.IsEqual(output, 'Chlorophyll-A'),
        imagesS3chlA,
        imagesS3
        )
    );
    animation.animate(images, {label: 'label', maxFrames: 100, compact: true, hidePlay: true, width: '400px'}); 
    
    //image collection for chart
    var imgCol = ee.Algorithms.If(
      ee.Algorithms.IsEqual(satellite, 'Sentinel-2'),
      S2filter,
      S3filter
      );
    var imgCol = ee.ImageCollection(imgCol);
    var reg = ee.Algorithms.If(
      ee.Algorithms.IsEqual(satellite, 'Sentinel-2'),
      SA,
      SAinner
      );
    var reg = ee.FeatureCollection(reg);
    
    //chl A chart for display
    var chlAChart = ui.Chart.image.doySeriesByYear({
      imageCollection: imgCol, 
      bandName: 'chlA',
      region: reg,
      regionReducer: ee.Reducer.mean(), 
      scale: 90
    });
    chlAChart.setOptions(chartOptions);
    var chartPanel = ui.Panel({widgets: [chlAChart],
      style: {
      stretch: 'horizontal', 
      backgroundColor: '#00000000',
      position: 'top-right',
      }
    });
    
    Map.add(chartPanel);
    Map.setOptions('HYBRID');
    Map.centerObject(SA, 12);
    Map.add(Panel);
    Map.add(legend);
    
    ///-----------------------------------------------------------------------------------------------------
    //re add About and intructions panel (for some reason this has to be all re added?)
    //------------------------------------------------------------------------------------------------------
    var CONTROL_PANEL_WIDTH = '280px';
    var CONTROL_PANEL_WIDTH_HIDE = '165px';
    var infoFont = {fontSize: '11px', color: '#505050'};
    var titleFont = {fontSize: '13px', fontWeight: 'bold', color: '#505050'};
    
    var instrucButton = ui.Button(
      {label: 'Instructions ❯', style: {margin: '0px 4px 0px 0px'}});
    
    var instrucTitle = ui.Label(
      'Instructions',
      titleFont);
    var instrucLabel1 = ui.Label(
      '1. In the Selection Panel, enter date range (2017-present day)',
      infoFont);
    var instrucLabel2 = ui.Label(
      '2. Map output will control whether users want to view chlorophyll A predictions or true colour RGB images',
      infoFont);
    var instrucLabel3 = ui.Label(
      '3. The Satellite tab allows users to select Sentinel-2 or Sentinel-3 data. ' + 
      'Sentinel-2 data is 20 m resolution but less frequent images while Sentinel-3 is 300 m resolution but more frequent.',
      infoFont);
    var instrucLabel4 = ui.Label(
      '4. Press “GO” to generate the images',
      infoFont);
    var instrucLabel5 = ui.Label(
      '5. Once the images have loaded (the progress will be shown in the top right layers tab), you can ' + 
      'click the animation bar and use arrow keys to scroll through the dates',
      infoFont);
    var instrucLabel6 = ui.Label(
      '6. A time-series chart will pop up for the desired date range.' + 
      'To download this data click the box at the top right of the chart and then select export to csv or png.',
      infoFont);
    
    var aboutTitle = ui.Label(
      'About',
      titleFont);
    var aboutButton = ui.Button(
      {label: 'About ❯', style: {margin: '0px 4px 0px 0px'}});
    var aboutLabel = ui.Label(
      'This web tool helps users visualize and download Sentinel-2 and Sentinel-3 derived chlorophyll A ' + 
      'predictions for Pigeon Lake, Alberta. It is a joint project between: the Alberta Biodiversity Monitoring ' +
      'Institute, University of Alberta, Alberta Lake Management Society, Alberta Environment and Parks, and ' +
      'Pigeon Lake Watershed Association. Chlorophyll A predictions are calibrated with fields samples ' +
      'obtained in spring/summer 2020. Contact edelance@ualberta.ca for questions.',
      infoFont);
      
    var controlPanel = ui.Panel({
      style: {position: 'bottom-right', width: CONTROL_PANEL_WIDTH_HIDE,
        maxHeight: '90%'
      }});
    
    var buttonPanel = ui.Panel(
      [instrucButton, aboutButton],
      ui.Panel.Layout.Flow('horizontal'),
      {stretch: 'horizontal', margin: '0px 0px 0px 0px'});
      
    var infoElements = ui.Panel(
      {style: {shown: false, margin: '0px -8px 0px -8px'}});
      
    var aboutElements = ui.Panel(
      {style: {shown: false, margin: '0px -8px 0px -8px'}});
      
    var infoShow = false;
    var controlShow = false;
    function infoButtonHandler() {
      if(infoShow) {
        infoShow = false;
        infoElements.style().set('shown', false);
        instrucButton.setLabel('Instructions ❯');
      } else {
        infoShow = true;
        infoElements.style().set('shown', true);
        instrucButton.setLabel('Instructions ❮');
      }
      if(infoShow | controlShow) {
        controlPanel.style().set('width', CONTROL_PANEL_WIDTH);
      } else {
        controlPanel.style().set('width', CONTROL_PANEL_WIDTH_HIDE);
      }
    }
    
    function aboutButtonHandler() {
      if(controlShow) {
        controlShow = false;
        aboutElements.style().set('shown', false);
        aboutButton.setLabel('About ❯');
      } else {
        controlShow = true;
        aboutElements.style().set('shown', true);
        aboutButton.setLabel('About ❮');
      }
      if(infoShow | controlShow) {
        controlPanel.style().set('width', CONTROL_PANEL_WIDTH);
      } else {
        controlPanel.style().set('width', CONTROL_PANEL_WIDTH_HIDE);
      }
    }
    
    infoElements.add(instrucTitle);
    infoElements.add(instrucLabel1);
    infoElements.add(instrucLabel2);
    infoElements.add(instrucLabel3);
    infoElements.add(instrucLabel4);
    infoElements.add(instrucLabel5);
    infoElements.add(instrucLabel6);
    
    aboutElements.add(aboutTitle);
    aboutElements.add(aboutLabel);
    
    
    
    controlPanel.add(buttonPanel);
    controlPanel.add(infoElements);
    controlPanel.add(aboutElements);
    
    
    Map.add(controlPanel);
    instrucButton.onClick(infoButtonHandler);
    aboutButton.onClick(aboutButtonHandler);
    //-----------------------------------------------------------------------------------------------------
    //
    //-----------------------------------------------------------------------------------------------------
    
  },
  style: {
    backgroundColor: C2
  }
});
//----------------------------------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------
//Add selection panel
//-----------------------------------------------------------------------------------------------------
var Panel = ui.Panel({
  widgets: [panelTitle,titleStartDate, textStartDate, titleEndDate, textEndDate, selectBox, selectBoxSensor, button],
  style: {
    position: 'top-left',
    backgroundColor: C1
  }
});
Map.add(Panel);
//-----------------------------------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------------------------------




