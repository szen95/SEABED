var classesNumber = 6,
    cellSize = 24;

//#########################################################
function heatmap_display(url, heatmapId, paletteName) {

    //==================================================
    // References
    // http://bl.ocks.org/Soylent/bbff6cc507dca2f48792
    // http://bost.ocks.org/mike/selection/
    // http://bost.ocks.org/mike/join/
    // http://stackoverflow.com/questions/9481497/understanding-how-d3-js-binds-data-to-nodes
    // http://bost.ocks.org/mike/miserables/
    // http://bl.ocks.org/ianyfchang/8119685

    //==================================================
    var tooltip = d3.select(heatmapId)
        .append("div")
        .style("position", "absolute")
        .style("visibility", "visible");

    //==================================================
    // http://bl.ocks.org/mbostock/3680958
    function zoom() {
        svg.attr("transform", "translate(" + d3.event.translate + ")scale(" + d3.event.scale + ")");
    }

    // define the zoomListener which calls the zoom function on the "zoom" event constrained within the scaleExtents
    var zoomListener = d3.behavior.zoom().scaleExtent([0.1, 3]).on("zoom", zoom);

    //==================================================
    var viewerWidth = $(document).width();
    var viewerHeight = $(document).height();
    var viewerPosTop = 300;
    var legendPosTop = 80;
    var viewerPosLeft = 350;

    var legendElementWidth = cellSize * 4;

    var colors = ['#444444', '#f96881', '#fdff26', '#7fbf7f', '#0066ff'];


    // http://bl.ocks.org/mbostock/3680999
    var svg;

    //==================================================
    d3.json(url, function (error, data) {


        var arr = data.data;
        var row_number = arr.length;
        var col_number = arr[0].length;

        var colorScale = d3.scale.quantize()
            .domain([0.0, 1.0, 2.0, 3.0, 4.0])
            .range(colors);

        svg = d3.select(heatmapId).append("svg")
            .attr("width", viewerWidth)
            .attr("height", viewerHeight)
            .call(zoomListener)
            .append("g")
            .attr("transform", "translate(" + viewerPosLeft + "," + viewerPosTop + ")");

        svg.append('defs')
            .append('pattern')
            .attr('id', 'diagonalHatch')
            .attr('patternUnits', 'userSpaceOnUse')
            .attr('width', 4)
            .attr('height', 4)
            .append('path')
            .attr('d', 'M-1,1 l2,-2 M0,4 l4,-4 M3,5 l2,-2')
            .attr('stroke', '#000000')
            .attr('stroke-width', 1);

        var rowSortOrder = false;
        var colSortOrder = false;

        var rowLabels = svg.append("g")
            .attr("class", "rowLabels")
            .selectAll(".rowLabel")
            .data(data.index)
            .enter().append("text")
            .text(function (d) {
                return d.count > 1 ? d.join("/") : d;
            })
            .attr("x", 0)
            .attr("y", function (d, i) {
                return (i * cellSize);
            })
            .style("text-anchor", "end")
            .attr("transform", function (d, i) {
                return "translate(-3," + cellSize / 1.5 + ")";
            })
            .attr("class", "rowLabel mono")
            .attr("id", function (d, i) {
                return "rowLabel_" + i;
            })
            .on('mouseover', function (d, i) {
                d3.select('#rowLabel_' + i).classed("hover", true);
            })
            .on('mouseout', function (d, i) {
                d3.select('#rowLabel_' + i).classed("hover", false);
            })
            .on("click", function (d, i) {
                rowSortOrder = !rowSortOrder;
                sortByValues("r", i, rowSortOrder);
                d3.select("#order").property("selectedIndex", 0);
            });

        var colLabels = svg.append("g")
            .attr("class", "colLabels")
            .selectAll(".colLabel")
            .data(data.columns)
            .enter().append("text")
            .text(function (d) {
                // d.shift();
                return d.count > 1 ? join("/") : d;
            })
            .attr("x", 0)
            .attr("y", function (d, i) {
                return (i * cellSize);
            })
            .style("text-anchor", "left")
            .attr("transform", function (d, i) {
                return "translate(" + cellSize / 2 + ", -3) rotate(-90) rotate(45, 0, " + (i * cellSize) + ")";
            })
            .attr("class", "colLabel mono")
            .attr("id", function (d, i) {
                return "colLabel_" + i;
            })
            .on('mouseover', function (d, i) {
                d3.select('#colLabel_' + i).classed("hover", true);
            })
            .on('mouseout', function (d, i) {
                d3.select('#colLabel_' + i).classed("hover", false);
            })
            .on("click", function (d, i) {
                colSortOrder = !colSortOrder;
                sortByValues("c", i, colSortOrder);
                d3.select("#order").property("selectedIndex", 0);
            });

        var row = svg.selectAll(".row")
            .data(data.data)
            .enter().append("g")
            .attr("id", function (d) {
                return d.idx;
            })
            .attr("class", "row");

        var j = 0;
        var heatMap = row.selectAll(".cell")
            .data(function (d) {
                j++;
                return d;
            })
            .enter().append("svg:rect")
            .attr("x", function (d, i) {
                return i * cellSize;
            })
            .attr("y", function (d, i, j) {
                return j * cellSize;
            })
            .attr("rx", 4)
            .attr("ry", 4)
            .attr("class", function (d, i, j) {
                return "cell bordered cr" + j + " cc" + i;
            })
            .attr("row", function (d, i, j) {
                return j;
            })
            .attr("col", function (d, i, j) {
                return i;
            })
            .attr("width", cellSize)
            .attr("height", cellSize)
            .style("fill", function (d) {
                if (d != null) return colorScale(d);
                else return "url(#diagonalHatch)";
            })
            .on('mouseover', function (d, i, j) {
                d3.select('#colLabel_' + i).classed("hover", true);
                d3.select('#rowLabel_' + j).classed("hover", true);

                

                var data_cell = d;
                
                if (data_cell == -1) { 
                    tooltip.html('<div class="heatmap_tooltip">' + 'No plot' + '</div>');
                    tooltip.style("visibility", "visible");
                    
                // } else if (data_cell == 0) {
                //     tooltip.html('<div class="heatmap_tooltip">' + 'Check for plot' + '</div>');
                //     tooltip.style("visibility", "visible");
                // } else {
                } else {
                    tooltip.html('<div class="heatmap_tooltip">' + 'Plot available' + '</div>');
                    tooltip.style("visibility", "visible");
                } 

                    
            })
            .on('mouseout', function (d, i, j) {
                d3.select('#colLabel_' + i).classed("hover", false);
                d3.select('#rowLabel_' + j).classed("hover", false);
                
                tooltip.style("visibility", "hidden");

            })
            .on("mousemove", function (d, i) {
                tooltip.style("top", (d3.event.pageY - 55) + "px").style("left", (d3.event.pageX - 60) + "px");
                
            })
            .on('click', function (d, i, j) {
                
                d3.select('#colLabel_' + i);
                d3.select('#rowLabel_' + j);

                var cols = labels.columns[i];
                var rows = labels.index[j];
                
                rows = rows.replace(/ *\([^)]*\) */g, "");
                cols = cols.replace(/ *\([^)]*\) */g, "");

                final_rows = rows.replace(/\(|\)/g, "");
                final_cols = cols.replace(/\(|\)/g, "");

                console.log(final_cols);

                var image_file = "website_s1/img/POET_Visualisation_" + final_rows + "_" + final_cols + ".pdf";
                openIMG(image_file);
                
            });

        function openIMG(d) {
            window.open(d);
            // document.getElementById('box').src = d;
            // console.log(d);
            // $("#modal01").show();
            // $("#box").show();
            

            // $("#box").click(function (e) {
            //     $("#modal01").hide();
            //     $("#box").hide();
            //   });
            
        }
        
        var legendElementText = [{ "label": "No differential response (i)", "value": 0 },
        { "label": "Sensitive to both MAPK and PI3K-AKT pathway inhibitors (ii)", "value": 1 },
        { "label": "Preferential MAPK pathway sensitivity (iii)", "value": 2 },
        { "label": "Preferential PI3K-AKT pathway sensitivity (iv)", "value": 3 },
        { "label": "Sensitive to either a MAPK pathway inhibitor or PI3K-AKT pathway inhibitor and vice versa (Divergent response, v)", "value": 4 }];


        var logo = svg.append("svg:image")
            .attr("class", "logo")
            .attr("xlink:href", "website_s1/img/SEABED_logo.png")
            .attr("x", legendElementWidth * 4)
            .attr("y", legendPosTop / 15)
            .attr("width", "20%")
            .attr("height", "20%");

        var legend = svg.append("g")
            .attr("class", "legend")
            .attr("transform", "translate(0,-370)")
            .selectAll(".legendElement")
            .data(legendElementText)
            .enter().append("g")
            .attr("class", "legendElement");


        legend.append("svg:rect")
            .attr("x", function (d, i) {
                return legendElementWidth * i;
            })
            .attr("y", legendPosTop)
            .attr("class", "cellLegend bordered")
            .attr("width", legendElementWidth)
            .attr("height", cellSize / 2)
            .style("fill", function (d, i) {
                return colors[i];
            });


        legend.append("text")
            .attr("class", "mono legendElement")
            .text(function (d) {
                return d.label;
            })
            .attr("x", function (d, i) {
                return legendElementWidth * i;
            })
            .attr("y", legendPosTop + cellSize)
            .call(wrap, 90);
        
        var akt_drug = svg.append("text")
            .attr("class", "drug_class")
            .text("AKT/PI3K inhibitors")
            .attr("x", function (d, i) {
                return i * cellSize;
            })
            .attr("y", -30)
            .attr("transform", function (d, i) {
                return "translate(" + cellSize / 2 + ", -3) rotate(-90) rotate(45, 0, " + (i * cellSize) + ")";
            })
        
        var mapk_drug = svg.append("text")
            .attr("class", "drug_class")
            .text("MAPK inhibitors")
            .attr("x", -150)
            .attr("y", 0)

        var fig_legend = svg.append("svg:image")
            .attr("class", "fig_legend")
            .attr("xlink:href", "website_s1/img/figure3_legend.png")
            .attr("x", legendElementWidth * 4)
            .attr("y", legendPosTop / 15)
            .attr("width", "40%")
            .attr("height", "40%");
        
        var fig_legend_texg = svg.append("text")
            .attr("class", "drug_class")
            .text("Legend")
            .attr("x", legendElementWidth * 5)
            .attr("y", legendPosTop / 10)

        
        function wrap(text, width) {
            text.each(function () {
                var text = d3.select(this),
                    words = text.text().split(/\s+/).reverse(),
                    word,
                    line = [],
                    lineNumber = 0,
                    lineHeight = 1.2, // ems
                    x = text.attr("x"),
                    y = text.attr("y"),
                    dy = text.attr("dy") ? text.attr("dy") : 0;
                tspan = text.text(null).append("tspan").attr("x", x).attr("y", y).attr("dy", dy + "em");
                while (word = words.pop()) {
                    line.push(word);
                    tspan.text(line.join(" "));
                    if (tspan.node().getComputedTextLength() > width) {
                        line.pop();
                        tspan.text(line.join(" "));
                        line = [word];
                        tspan = text.append("tspan").attr("x", x).attr("y", y).attr("dy", ++lineNumber * lineHeight + dy + "em").text(word);
                    }
                }
            });
        }

        //==================================================
        // Change ordering of cells
        function sortByValues(rORc, i, sortOrder) {
            var t = svg.transition().duration(1000);
            var values = [];
            var sorted;
            d3.selectAll(".c" + rORc + i)
                .filter(function (d) {
                    if (d != null) values.push(d);
                    else values.push(-999); // to handle NaN
                });
            	
            if (rORc == "r") { // sort on cols
                sorted = d3.range(col_number).sort(function (a, b) {
                    if (sortOrder) {
                        return values[b] - values[a];
                    } else {
                        return values[a] - values[b];
                    }
                });
                t.selectAll(".cell")
                    .attr("x", function (d) {
                        var col = parseInt(d3.select(this).attr("col"));
                        return sorted.indexOf(col) * cellSize;
                    });
                t.selectAll(".colLabel")
                    .attr("y", function (d, i) {
                        return sorted.indexOf(i) * cellSize;
                    })
                    .attr("transform", function (d, i) {
                        return "translate(" + cellSize / 2 + ", -3) rotate(-90) rotate(45, 0, " + (sorted.indexOf(i) * cellSize) + ")";
                    });
            } else { // sort on rows
                sorted = d3.range(row_number).sort(function (a, b) {
                    if (sortOrder) {
                        return values[b] - values[a];
                    } else {
                        return values[a] - values[b];
                    }
                });
                t.selectAll(".cell")
                    .attr("y", function (d) {
                        var row = parseInt(d3.select(this).attr("row"));
                        return sorted.indexOf(row) * cellSize;
                    });
                t.selectAll(".rowLabel")
                    .attr("y", function (d, i) {
                        return sorted.indexOf(i) * cellSize;
                    })
                    .attr("transform", function (d, i) {
                        return "translate(-3," + cellSize / 1.5 + ")";
                    });
            }
        }

        //==================================================
        d3.select("#order").on("change", function () {
            var newOrder = d3.select("#order").property("value");
            changeOrder(newOrder, heatmapId);
        });

        //==================================================
        d3.select("#palette")
            .on("keyup", function () {
                var newPalette = d3.select("#palette").property("value");
                if (newPalette != null)						// when interfaced with jQwidget, the ComboBox handles keyup event but value is then not available ?
                    changePalette(newPalette, heatmapId);
            })
            .on("change", function () {
                var newPalette = d3.select("#palette").property("value");
                changePalette(newPalette, heatmapId);
            });
    });

    //==================================================
}

//#########################################################
function changeOrder(newOrder, heatmapId) {
    var svg = d3.select(heatmapId);
    var t = svg.transition().duration(1000);
    if (newOrder == "sortinit_col") { // initial sort on cols (alphabetically if produced like this)
        t.selectAll(".cell")
            .attr("x", function (d) {
                var col = parseInt(d3.select(this).attr("col"));
                return col * cellSize;
            });
        t.selectAll(".colLabel")
            .attr("y", function (d, i) {
                return i * cellSize;
            })
            .attr("transform", function (d, i) {
                return "translate(" + cellSize / 2 + ", -3) rotate(-90) rotate(45, 0, " + (i * cellSize) + ")";
            });
    } else if (newOrder == "sortinit_row") { // initial sort on rows (alphabetically if produced like this)
        t.selectAll(".cell")
            .attr("y", function (d) {
                var row = parseInt(d3.select(this).attr("row"));
                return row * cellSize;
            });
        t.selectAll(".rowLabel")
            .attr("y", function (d, i) {
                return i * cellSize;
            })
            .attr("transform", function (d, i) {
                return "translate(-3," + cellSize / 1.5 + ")";
            });
    } else if (newOrder == "sortinit_col_row") { // initial sort on rows and cols (alphabetically if produced like this)
        t.selectAll(".cell")
            .attr("x", function (d) {
                var col = parseInt(d3.select(this).attr("col"));
                return col * cellSize;
            })
            .attr("y", function (d) {
                var row = parseInt(d3.select(this).attr("row"));
                return row * cellSize;
            });
        t.selectAll(".colLabel")
            .attr("y", function (d, i) {
                return i * cellSize;
            })
            .attr("transform", function (d, i) {
                return "translate(" + cellSize / 2 + ", -3) rotate(-90) rotate(45, 0, " + (i * cellSize) + ")";
            });
        t.selectAll(".rowLabel")
            .attr("y", function (d, i) {
                return i * cellSize;
            })
            .attr("transform", function (d, i) {
                return "translate(-3," + cellSize / 1.5 + ")";
            });
    }
}

//#########################################################
function changePalette(paletteName, heatmapId) {
    var colors = colorbrewer[paletteName][classesNumber];
    var colorScale = d3.scale.quantize()
        .domain([0.0, 1.0])
        .range(colors);
    var svg = d3.select(heatmapId);
    var t = svg.transition().duration(500);
    t.selectAll(".cell")
        .style("fill", function (d) {
            if (d != null) return colorScale(d);
            else return "url(#diagonalHatch)";
        })
    t.selectAll(".cellLegend")
        .style("fill", function (d, i) {
            return colors[i];
        });
}
