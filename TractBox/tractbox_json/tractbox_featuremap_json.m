function rootjson = tractbox_featuremap_json(options,type,sourcesize,dotimport)
rootjson.description = 'Feature map as described in Lambert et al., 2017 Defining thalamic nuclei and topographic connectivity gradients in vivo';
rootjson.tractographyinput = dotimport;
rootjson.type=type;
rootjson.sourcesize = sourcesize;
rootjson.threshold = options.featuremap.thr;
end