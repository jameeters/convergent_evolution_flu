"""
Get a list of all avian hostnames
"""
avians = set()
with open('data/all_hosts.txt', 'r') as f:
    hosts = [l.strip() for l in f.readlines()]

# Avian keywords
avian_keywords = [
    'avian', 'bird', 'chicken', 'duck', 'goose', 'turkey', 'hawk', 'eagle', 'owl', 
    'swan', 'teal', 'gull', 'wigeon', 'pintail', 'crane', 'pelican', 'vulture', 
    'grebe', 'scoter', 'merganser', 'tern', 'sandpiper', 'heron', 'egret', 
    'quail', 'pheasant', 'dove', 'pigeon', 'crow', 'raven', 'magpie', 'jay', 
    'starling', 'sparrow', 'robin', 'thrush', 'blackbird', 'grackle', 'falcon', 
    'osprey', 'kite', 'harrier', 'coot', 'stilt', 'shoveler', 'gadwall', 'mallard', 
    'canvasback', 'redhead', 'scaup', 'brant', 'eider', 'goldeneye', 'bufflehead', 
    'merlin', 'kestrel', 'cormorant', 'gannet', 'guinea', 'peafowl', 'emu', 
    'flamingo', 'ibis', 'plover', 'turnstone', 'dunlin', 'willet', 'dowitcher', 
    'yellowlegs', 'phalarope', 'murre', 'puffin', 'auk', 'loon', 'ostrich', 
    'rhea', 'kiwi', 'gallus', 'meleagris', 'anas', 'anser', 'aythya', 'branta', 
    'buteo', 'corvus', 'cygnus', 'falco', 'larus', 'numididae', 'phasianus', 
    'spatula', 'strigiformes', 'strix', 'tyto', 'bubo', 'accipiter', 'haliaeetus', 
    'cathartes', 'pandion', 'gavia', 'pelecanus', 'ardea', 'egretta', 'cairina', 
    'aix', 'bucephala', 'callipepla', 'calidris', 'albatross', 'petrel', 'condor',
    'raptor', 'fowl', 'poultry', 'parakeet', 'parrot', 'cockatoo', 'finch', 'bunting',
    'warbler', 'cardinal', 'junco', 'chickadee', 'nuthatch', 'wren', 'mockingbird',
    'catbird', 'swallow', 'martin', 'swift', 'hummingbird', 'kingfisher', 'woodpecker',
    'flycatcher', 'shrike', 'vireo', 'lark', 'bluebird', 'waxwing', 'towhee', 'oriole',
    'goldfinch', 'siskin', 'redpoll', 'crossbill', 'tundra', 'trumpeter', 'shearwater',
    'booby', 'anhinga', 'bittern', 'stork', 'limpkin', 'rail', 'gallinule', 'oystercatcher',
    'avocet', 'lapwing', 'godwit', 'whimbrel', 'curlew', 'snipe', 'woodcock', 'skua', 'noddy',
    'skimmer', 'alcid', 'cuckoo', 'roadrunner', 'trogon', 'sapsucker', 'flicker', 'phoebe',
    'kingbird', 'gnatcatcher', 'dipper', 'kinglet', 'solitaire', 'thrasher', 'meadowlark',
    'cowbird', 'weaver', 'whydah', 'indigobird', 'black scoter', 'surf scoter', 'white-winged scoter',
    'bufflehead', 'common goldeneye', 'barrow\'s goldeneye', 'ruddy duck', 'wood duck', 
    'fulvous whistling-duck', 'black-bellied whistling-duck', 'northern bobwhite', 'wild turkey',
    'ruffed grouse', 'spruce grouse', 'willow ptarmigan', 'rock ptarmigan', 'white-tailed ptarmigan',
    'dusky grouse', 'sooty grouse', 'sharp-tailed grouse', 'greater prairie-chicken', 
    'lesser prairie-chicken', 'ring-necked pheasant', 'gray partridge', 'chukar',
    'lophodytes', 'accipitridae', 'thalasseus', 'phasianidae', 'mergus', 'mareca', 'arinae', 
    'laridae', 'antigone', 'sibirionetta', 'peacock', 'peahen', 'cago', 'hosp', 'pefa',
    'caracara', 'sanderling', 'melegris', 'corvidae', 'aves', 'pelecanidae', 'phalacrocoracidae',
    'partridge', 'peregrine', 'kittiwake', 'gadwell', 'widgeon', 'herrier'
]

for host in hosts:
    if any(k in host for k in avian_keywords):
        avians.add(host)

print(avians)