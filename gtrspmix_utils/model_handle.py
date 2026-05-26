
import re

def parse_model(nex): # read nexus file
    l_freq, l_fmix, l_exch, l_mix = [],[],[],[]
    text = nex.read_text()
    text = re.sub('\n{2,}', '\n', text)
    text = re.sub('\n +', '\n', text)
    text = re.sub('; +', ';', text)
    text = re.sub(' {2,}', ' ', text)
    models = text.split(';\n')[1:-2]
    for model in models:
        model_type = model.split(' ')[0]
        model_name = model.split(' ')[1].split('=')[0].replace(' ', '')
        model_values = model.rstrip().split('=')[1].lstrip(' ')
        if model_type == 'frequency' and not 'FMIX' in model_values:
            l_freq.append([model_name, model_values])
        elif model_type == 'frequency' and 'FMIX' in model_values:
            l_fmix.append([model_name, model_values])
        elif model_type == 'model' and not 'MIX' in model_values:
            l_exch.append([model_name, model_values])
        elif model_type == 'model' and 'MIX' in model_values:
            l_mix.append([model_name, model_values])

    return l_freq, l_fmix, l_exch, l_mix


def write_model(outfile, l_freq, l_fmix, l_exch, l_mix): # write nexus file
    header = '#nexus\nbegin models;\n'
    footer = '\nend;\n'
    with outfile.open(mode='w') as fo:
        print(header, file=fo)
        for j in [l_freq, l_fmix]:
            if len(j) > 0:
                for i in j:
                    print(f'frequency {i[0]} = {i[1]} ;', file=fo)
                print('', file=fo)
        for j in [l_exch, l_mix]:
            if len(j) > 0:
                for i in j:
                    print(f'model {i[0]} = {i[1]} ;', file=fo)
                print('', file=fo)
        print(footer, file=fo)


def make_fmix(d_cluster, l_freq):
    if d_cluster == None:
        l_freq_comp = []
        for freq in l_freq:
            l_freq_comp.append(freq[0])
        model = f'FMIX{{{",".join(l_freq_comp)}}}'
        l_fmix = [[f'Freq{len(l_freq)}', model]]

    else:
        l_fmix = []
        for i in d_cluster:
            l_cluster_freq = []
            for j in d_cluster[i]:
                l_cluster_freq.append(l_freq[j-1][0])
            model = f'FMIX{{{",".join(l_cluster_freq)}}}'
            l_fmix.append([f'Fcluster{i}', model])
    return l_fmix


def init_model(args, d_cluster, l_freq, l_exch): # initialize model structure
    l_fmix = make_fmix(d_cluster, l_freq)

    l_all_mix = []
    for i in l_freq:
        l_all_mix.append(f'{args.m_gtr20}+F{i[0]}')
    model = f'MIX{{{",".join(l_all_mix)}}}'
    num_exch = len(d_cluster)
    l_mix = [[f'S{num_exch}F{len(l_freq)}Opt', model]]

    if l_exch == []: # create dummy l_exch
        for i in range(num_exch):
            l_exch.append([args.m_gtr20, ''])

    return l_fmix, l_exch, l_mix


def extract_LL_rate_modelw(iqtree_file): # extract params from iqtree file
    iqtree_file = iqtree_file.read_text()
    LL = iqtree_file.split('Log-likelihood of the tree: ', 1)[1].split(' ', 1)[0]

    rate_type = iqtree_file.split('Model of rate heterogeneity: ', 1)[1].split(' ', 1)[0]
    if rate_type == 'Gamma':
        rate_param = iqtree_file.split('Gamma shape alpha: ', 1)[1].split('\n', 1)[0]
    elif rate_type == 'FreeRate':
        rate_param = iqtree_file.split('Site proportion and rates:  ', 1)[1].split('\n', 1)[0]
        rate_param = rate_param.replace(') (', ',')[1:-1]

    components = iqtree_file.split('SUBSTITUTION PROCESS', 1)[1]
    components = components.split('Parameters\n', 1)[1]
    components = components.split('\n\n', 1)[0]

    l_mix_comp = []
    d_fmix_comp = {}
    for i in components.split('\n'):
        l = i.split()
        comp = f'{l[1]}:{float(l[2])}:{l[3]}'
        l_mix_comp.append(comp)
        d_fmix_comp[l[1].split('+F')[1]] = f':{float(l[2])}:{l[3]}'
    l_mix_w = 'MIX{' + ','.join(l_mix_comp) + '}'
    return float(LL), rate_param, l_mix_w, d_fmix_comp


def extract_GTR(args, d_cluster, l_freq, num_it): # extract GTRs from iqtree file
    new_l_exch = []
    for i in d_cluster:
        iqtree_file = args.outdir / f'GTR_opt_{num_it}' / f'cluster{i}' / 'input_train.fasta.iqtree'
        exch = iqtree_file.read_text().split('Linked substitution parameters (lower-diagonal):\n')[1].split('\n\nModel of rate heterogeneity:')[0]
        exch = exch + '\n' + ' 0.05'*20
        new_l_exch.append([f'S{len(d_cluster)}S{i}', exch])
    l_all_comp = []
    for i in range(len(l_freq)):
        for j in d_cluster:
            if i+1 in d_cluster[j]:
                l_all_comp.append(f'S{len(d_cluster)}S{j}+F{l_freq[i][0]}')
    all_comp = f'MIX{{{",".join(l_all_comp)}}}'
    new_l_mix = [[f'S{len(d_cluster)}F{len(l_freq)}Opt', all_comp]]
    return new_l_exch, new_l_mix


def extract_Profile(args, num_it, d_cluster, l_freq, l_exch): # FO
    rev = {x: k for k, lst in d_cluster.items() for x in lst}
    new_l_freq = []
    for i in range(len(l_freq)):
        iqtree_file = args.outdir / f'Profile_opt_{num_it}' / f'p{i+1}' / 'input_train.fasta.iqtree'
        text = iqtree_file.read_text().split('State frequencies: (estimated with maximum likelihood)\n\n')[1].split('\n\nModel of rate heterogeneity:')[0] # FO
        # text = iqtree_file.read_text().split('State frequencies: (empirical counts from alignment)\n\n')[1].split('\n\nModel of rate heterogeneity:')[0] # F
        l_profile = []
        for l in text.split('\n'):
            l_profile.append(l.split(' = ')[1])
        profile = ' '.join(l_profile)
        new_l_freq.append([f'F{len(l_freq)}pi{i+1}', profile])

    # l_fmix
    new_l_fmix = make_fmix(d_cluster, new_l_freq)

    # l_mix
    l_mix_comp = []
    for i in range(len(l_freq)):
        l_mix_comp.append(f'{l_exch[int(rev[i+1])-1][0]}+F{new_l_freq[i][0]}')
    model = f'MIX{{{",".join(l_mix_comp)}}}'
    num_exch = len(d_cluster)
    new_l_mix = [[f'S{num_exch}F{len(l_freq)}Opt', model]]

    return new_l_freq, new_l_fmix, new_l_mix


def get_Profile_F(profiles, d_cluster, l_freq, l_exch): # F
    rev = {x: k for k, lst in d_cluster.items() for x in lst}
    new_l_freq = []
    for m in range(profiles.shape[1]):
        profile = profiles[:, m]
        l_profile = [f"{freq:.10f}" for freq in profile]
        profile = ' '.join(l_profile)
        new_l_freq.append([f'F{len(l_freq)}pi{m+1}', profile])

    # l_fmix
    new_l_fmix = make_fmix(d_cluster, new_l_freq)

    # l_mix
    l_mix_comp = []
    for i in range(len(l_freq)):
        l_mix_comp.append(f'{l_exch[int(rev[i+1])-1][0]}+F{new_l_freq[i][0]}')
    model = f'MIX{{{",".join(l_mix_comp)}}}'
    num_exch = len(d_cluster)
    new_l_mix = [[f'S{num_exch}F{len(l_freq)}Opt', model]]

    return new_l_freq, new_l_fmix, new_l_mix
