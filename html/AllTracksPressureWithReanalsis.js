function makeplot(gdlid) {
  Plotly.d3.csv("https://raw.githubusercontent.com/Rafnuss/GeoPressureMAT/main/data/export/label/AllTracksPressureWithReanalysis_"+gdlid+".csv", function(data){ 

    traces = data.reduce(function (acc, d) {
      if (d.series == "presGL") {
        acc[0].x.push(d.date);
        acc[0].y.push(d.obs);
        if (d.isOutliar=="1"){
          acc[1].x.push(d.date);
          acc[1].y.push(d.obs);
        }
        
      } else {
        acc[2].x.push(d.date);
        acc[2].y.push(d.obs);
        //acc[1].outliar.push(d.isOutliar);
      }
      return acc
    }, [{
        name: 'Geolocator',
        x:[],
        y:[]
      }, 
      {
        name: 'Geolocator excluded',
        x:[],
        y:[],
        mode: 'markers',
        type: 'scatter'
      }, 
      {
        name: 'ERA5 best match',
        x:[],
        y:[]
      }]);

    
    Plotly.newPlot('plotlyID', traces,{
      xaxis: {
        autorange: true,
        rangeselector: {buttons: [
            {
              count: 1,
              label: '1d',
              step: 'day',
              stepmode: 'backward'
            },
            {
              count: 1,
              label: '6m',
              step: 'month',
              stepmode: 'backward'
            },
            {step: 'all'}
          ]},
        type: 'date'
      },
      yaxis: {
        autorange: true,
        type: 'linear'
      },
      legend: {
        x: 1,
        xanchor: 'right',
        y: 1,
        "orientation": "h"
      },
      hovermode: 'none',
    });

  });
};

makeplot("16AQ");

document.getElementById('track').onchange = function () {
  makeplot(this.value)
}