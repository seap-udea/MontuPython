# Montuheka parameters
montuheka_params = dict(
    angdist = [1,0,180],
    angspeed = [1,0,1],
    sundec = [1,-24,2*24] 
)

# Normalize montuheka function parameters
pars = []
i = 0
keys = dict()
for key,item in montuheka_params.items():
    pars += [item[0]]
    keys[i] = key
    i += 1
pars,norm = spy.unorm(pars)
for i in range(len(montuheka_params.keys())):
    montuheka_params[keys[i]][0] = pars[i]

# Montuheka function
def montuheka_function(series):
    mhf = sum([montuheka_params[prop][0]*\
         (series[prop]-montuheka_params[prop][1])/\
            montuheka_params[prop][2] for prop in montuheka_params.keys()])
    return mhf