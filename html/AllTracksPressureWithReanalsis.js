function makeplot(gdlid) {
  document.getElementById('loader').style.visibility = "visible";
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
        staId = acc.findIndex(x=>x.staid==d.sta_id)
        if (staId==-1){
          acc.push({
            name: "ERA5",
            staid: d.sta_id,
            x:[],
            y:[],
            mode: 'lines',
            showlegend: d.sta_id=="1"
          })
          staId = acc.findIndex(x=>x.staid==d.sta_id)
        } 
        acc[staId].x.push(d.date);
        acc[staId].y.push(d.obs);
      }
      return acc
    }, [{
        name: 'Geolocator',
        x:[],
        y:[],
        line: {
          color: 'grey',
          width: 3
        }
      }, 
      {
        name: 'Geolocator excluded',
        x:[],
        y:[],
        mode: 'markers',
        type: 'scatter',
        marker: {
          color:"black"
        }
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
      hovermode: false,
    });
    document.getElementById('loader').style.visibility = "hidden";
  });
};

makeplot("16AQ");

document.getElementById('track').onchange = function () {
  makeplot(this.value)
}