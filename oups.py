trouts = {"IHWH" : "IMG:3300020735",
"IHXT" : "IMG:3300020721",
"IHUI" : "IMG:3300020685"}

all_samples.update( {s : {'dataset' : 'troutbog', 'accession' : trouts[s]} for s in samples3})
