% Helper function with presets.
function colMap = customcolormap_preset(preset_name)

    switch(preset_name)
        case 'pasteljet'
            colMap = customcolormap([0 .25 .5 .75 1], {'#9d0142','#f66e45','#ffffbb','#65c0ae','#5e4f9f'});
        case 'red-white-blue'
            colMap = customcolormap(linspace(0,1,11), {'#68011d','#b5172f','#d75f4e','#f7a580','#fedbc9','#f5f9f3','#d5e2f0','#93c5dc','#4295c1','#2265ad','#062e61'});
        case 'white-red'
            colMap = customcolormap(linspace(0,1,6), {'#68011d', '#b5172f', '#d75f4e', '#f7a580', '#fedbc9', '#f5f9f3'});
        otherwise
            error('Unknown preset.');
    end
end