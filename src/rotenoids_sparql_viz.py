import pandas as pd
from pandas import json_normalize
from SPARQLWrapper import SPARQLWrapper, JSON
import pandas as pd
import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from tqdm import tqdm
import plotly.express as px
import tmap as tm
import numpy as np
from map4 import MAP4Calculator
import datatable as dt
from tqdm import tqdm
import umap
import matplotlib.pyplot as plt

uri_ref = 'https://enpkg.commons-lab.org/kg/'
uri_module = 'https://enpkg.commons-lab.org/module/'

sparql = SPARQLWrapper('http://FASIE-1439:7200/repositories/graph_sandbox')
sparql.setReturnFormat(JSON)

### TMAP OF MEMO VECTORS OF THE 1,600 SET


# Get samples metadata
sparql.setQuery("""
    PREFIX enpkg: <https://enpkg.commons-lab.org/kg/>
    PREFIX enpkgmodule: <https://enpkg.commons-lab.org/module/>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    PREFIX wdt: <http://www.wikidata.org/prop/direct/>
    
    SELECT ?extract ?taxon ?organe ?tcact ?tox
    
    WHERE
    { 
        ?extract rdf:type enpkg:LabExtract .
          ?extract enpkgmodule:has_bioassay_results ?biores .
          ?extract enpkgmodule:has_bioassay_results ?toxres .
            ?biores rdf:type enpkgmodule:Tcruzi10ugml .
            ?toxres rdf:type enpkgmodule:L610ugml .
              ?biores enpkgmodule:inhibition_percentage ?tcact .
              ?toxres enpkgmodule:inhibition_percentage ?tox .
         ?material enpkg:has_lab_process ?extract .
            ?material enpkgmodule:has_organe ?organe
            OPTIONAL {
              ?material enpkg:has_wd_id ?wd_id .
              service <https://query.wikidata.org/sparql> {
                  ?wd_id wdt:P225 ?taxon
                }
              }
    }
"""
)

results = sparql.queryAndConvert()
df = json_normalize(results['results']["bindings"])
df = df.stack().str.replace(uri_ref, "").unstack()
df.drop(list(df.filter(regex = 'type')), axis = 1, inplace = True)
df.columns = df.columns.str.replace('.value', '')
df.taxon = df.taxon.fillna('unresolved')
df.tcact = df.tcact.astype(float)
df.tox = df.tox.astype(float)

df['active'] = np.where(((df.tcact >= 80)) & (df.tox <= 50), 'Active & non-toxic', 'Inactive and/or toxic')

#Load MEMO analysis

def compute_tmap(table):
    
    dims = 1024
    enc = tm.Minhash(len(table.columns), 42, dims)

    i = 0
    fps = []
    for _, row in tqdm(table.iterrows()):
        fps.append(tm.VectorFloat(list(row))) #VectorFloat or VectorUint
        i+=1

    lf = tm.LSHForest(dims * 2, 128, store=True, weighted=True)
    lf.batch_add(enc.batch_from_weight_array(fps, method="I2CWS")) #int
    lf.index()

    CFG_TMAP = tm.LayoutConfiguration()
    CFG_TMAP.k = 10
    CFG_TMAP.kc = 10

    x, y, s, t, _ = tm.layout_from_lsh_forest(lf, CFG_TMAP)
    return x, y, s, t

memo_matrix = dt.fread('../data/memo_both_blanks_occ_0.gz').to_pandas()
memo_matrix.set_index('filename', inplace=True)
memo_matrix = memo_matrix[memo_matrix.index.isin(df.extract)]
memo_matrix = memo_matrix.reindex(df.extract)

x, y, s, t = compute_tmap(memo_matrix)

#del(memo_matrix)

df['tmap_x'] = x
df['tmap_y'] = y

px.scatter(df, x = 'tmap_x', y = 'tmap_y', color = 'active', hover_name='taxon', hover_data= ['organe'])

reducer = umap.UMAP(metric='braycurtis')
embedding = reducer.fit_transform(memo_matrix.values)
df['umap_x'] = embedding[:, 0]
df['umap_y'] = embedding[:, 1]
px.scatter(df, x = 'umap_x', y = 'umap_y', color = 'active',  hover_name='taxon', hover_data= ['organe'])

categories = ['active']
dic_cat= {}
for cat in categories:
    if cat == 'active':
        categorical = True
        colorsIdx = {
        'Active & non-toxic': '#06d6a0',
        'Inactive and/or toxic': '#e9d8a6',
        }
        title = '<i>T. cruzi<\i> activity'

        
    dic_cat[cat] = {}
    dic_cat[cat]['categorical'] = categorical
    dic_cat[cat]['colorsIdx'] = colorsIdx
    dic_cat[cat]['title'] = title
  
  
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=False, constrained_layout=True)
fig.set_size_inches(30, 11)
fig.set_dpi(300)
fig.suptitle("TMAP and UMAP vizualisations of samples' MEMO vectors colored according to their anti-$\it{T. cruzi}$ activity", fontsize=30)

df['size'] = np.where((df['active'] == 'Active & non-toxic'), 120, 40)  

cat = 'active'
for ax, coord in zip([ax1, ax2, ax3], [('tmap_x', 'tmap_y'), None, ('umap_x', 'umap_y')]):
  if ax == ax2:
    ax.axis('off')
  else:
    colors = list(map(dic_cat[cat]['colorsIdx'].get, list(df[cat])))
    for activity in dic_cat[cat]['colorsIdx'].keys():
      df_sel = df[df[cat] == activity]
      color = dic_cat[cat]['colorsIdx'][activity]

      if activity == 'Active & non-toxic':
        zorder = 3
      else:
        zorder=2
      
      scatter = ax.scatter(df_sel[coord[0]], df_sel[coord[1]], s=df_sel['size'], c=color, label=activity,
                          alpha = 0.8, edgecolors='grey',zorder=zorder)
      legend = ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.3),
          ncol=1, fancybox=False, fontsize=22)
      if ax == ax3:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        # Only show ticks on the left and bottom spines
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
        ax.tick_params(axis='x', labelsize=20)
        ax.tick_params(axis='y', labelsize=20)
        ax.set_xlabel('X', fontsize=30)
        ax.set_ylabel('Y', fontsize=30)
        
      else:
        ax.axis('off')
      if ax == ax1:
        for i in range(len(s)):
          ax.plot(
              [x[s[i]], x[t[i]]],
              [y[s[i]], y[t[i]]],
              "k-",
              alpha=0.5,
              zorder=1,    
          )
            
plt.savefig('../output/rotenoids/memo_both_plot_activity.png', bbox_inches='tight')

### TMAP OF ANNOTATIONS OF SELECTED SAMPLES

selection = 'VGF150_H08|VGF146_F09|VGF146_E08|VGF150_A09|VGF155_D05|VGF139_H07'

# Return annotations in selected samples, with their count 
sparql.setQuery("""
    PREFIX enpkg: <https://enpkg.commons-lab.org/kg/>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    
    SELECT DISTINCT ?ik_2d ?smiles ?np_class ?np_superclass ?np_pathway ?count ?count_in_group
    WHERE {
      ?ik_2d rdf:type ?InChIkey2D . 
        ?ik_2d enpkg:has_smiles ?smiles . 
        ?ik_2d enpkg:has_npc_class ?np_class .
        ?ik_2d enpkg:has_npc_superclass ?np_superclass .
        ?ik_2d enpkg:has_npc_pathway ?np_pathway .
          ?annotation enpkg:has_InChIkey2D ?ik_2d . 
              ?feature enpkg:has_sirius_annotation|enpkg:has_isdb_annotation ?annotation .
                ?feature_list enpkg:has_lcms_feature ?feature .
                  ?lcms enpkg:has_lcms_feature_list ?feature_list .
                    ?sample enpkg:has_LCMS ?lcms .
                    FILTER(regex(str(?sample), "%s"))
    {
      SELECT DISTINCT ?ik_2d (COUNT(DISTINCT(?sample)) AS ?count_in_group)
      WHERE {   
        {?annotation rdf:type enpkg:SiriusStructureAnnotation .
            ?annotation enpkg:has_cosmic_score ?cosmic .
            FILTER(?cosmic >= 0.4)}
        UNION
        {?annotation rdf:type enpkg:IsdbAnnotation .
          ?annotation enpkg:has_isdb_taxo_score ?taxo_score .
          FILTER(?taxo_score >= 6)} 
        ?annotation enpkg:has_InChIkey2D ?ik_2d . 
          ?feature enpkg:has_sirius_annotation|enpkg:has_isdb_annotation ?annotation .
            ?feature_list enpkg:has_lcms_feature ?feature .
              ?lcms enpkg:has_lcms_feature_list ?feature_list .
                  ?sample enpkg:has_LCMS ?lcms .
                  FILTER(regex(str(?sample), "%s"))
      } GROUP BY ?ik_2d ORDER BY DESC(?count_in_group)
    }
      {
        SELECT DISTINCT ?ik_2d (COUNT(DISTINCT(?sample)) AS ?count)
        WHERE {   
          {?annotation rdf:type enpkg:SiriusStructureAnnotation .
            ?annotation enpkg:has_cosmic_score ?cosmic .
            FILTER(?cosmic >= 0.4)}
          UNION
          {?annotation rdf:type enpkg:IsdbAnnotation .
            ?annotation enpkg:has_isdb_taxo_score ?taxo_score .
            FILTER(?taxo_score >= 6)} 
          ?annotation enpkg:has_InChIkey2D ?ik_2d . 
            ?feature enpkg:has_sirius_annotation|enpkg:has_isdb_annotation ?annotation .
              ?feature_list enpkg:has_lcms_feature ?feature .
                 ?lcms enpkg:has_lcms_feature_list ?feature_list .
                   ?sample enpkg:has_LCMS ?lcms .
        } GROUP BY ?ik_2d  
      }
      #FILTER((?count <= 80))
    }
""" % (selection, selection)
)

results = sparql.queryAndConvert()
df = json_normalize(results['results']["bindings"])
df = df.stack().str.replace(uri_ref, "").unstack()
df.drop(list(df.filter(regex = 'type')), axis = 1, inplace = True)
df.columns = df.columns.str.replace('.value', '')
df = df[df['smiles'] != 'unknown']

#Plotting t-SNE

def fp_list_from_smiles_list(smiles_list,n_bits=2048):
    fp_list = []
    for smiles in tqdm(smiles_list):
        mol = Chem.MolFromSmiles(smiles)
        fp_list.append(fp_as_array(mol,n_bits))
    return fp_list

def fp_as_array(mol,n_bits=2048):
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=n_bits)
    arr = np.zeros((1,), np.int)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr

fp_list = fp_list_from_smiles_list(df.smiles)
pca = PCA(n_components=2)
crds = pca.fit_transform(fp_list) 
crds_df = pd.DataFrame(crds,columns=["PC_1","PC_2"])
df["PC_1"] = list(crds_df["PC_1"])
df["PC_2"] = list(crds_df["PC_2"])

pca = PCA(n_components=50)  
crds_embedded = TSNE(n_components=2).fit_transform(crds)
tsne_df = pd.DataFrame(crds_embedded,columns=["X","Y"])
df["X"] = list(tsne_df["X"])
df["Y"] = list(tsne_df["Y"])
df['size'] = 1/df['count'].astype(int)
df['count_in_group']=df['count_in_group'].astype(int)

fig = px.scatter(df, x = 'X', y ='Y', color = 'count_in_group', hover_name='np_class', color_continuous_scale='Viridis')
fig.update_layout(template='simple_white', )
fig.show()


# Plot TMAP annotations

def compute_tmap(table):
    
    #enc = MHFPEncoder(1024)
    MAP4 = MAP4Calculator(dimensions=1024)

    lf = tm.LSHForest(1024, 64)
    
    i = 0
    fps = []
    for _, row in table.iterrows():
        if i != 0 and i % 10 == 0:
                print(100 * i / len(table))
        mol = AllChem.MolFromSmiles(row["smiles"])
        fps.append(mol)        
        i+=1

    fps = MAP4.calculate_many(fps)
    lf.batch_add(fps)
    lf.index()

    CFG_TMAP = tm.LayoutConfiguration()
    CFG_TMAP.k = 20
    CFG_TMAP.kc = 20

    x, y, s, t, _ = tm.layout_from_lsh_forest(lf, CFG_TMAP)
    return x, y, s, t
  
x, y, s, t = compute_tmap(df)

df['x'] = x
df['y'] = y
df['count_in_group']=df['count_in_group'].astype(int)
df['sel_class'] = np.where(((df['np_class'] == 'npc_Rotenoids')|(df['np_class'] == 'npc_Isoflavones')), df['np_class'], 'Other')
df['specificity'] = df['count_in_group'].astype(int)/df['count'].astype(int)
df['sel_pathway'] = np.where((
  (df['np_pathway'] == 'npc_Terpenoids')|
  (df['np_pathway'] == 'npc_Shikimates_and_Phenylpropanoids')|
  (df['np_pathway'] == 'npc_Fatty_acids')), df['np_pathway'], 'Other')

px.scatter(x = df['x'], y = df['y'], color = df['sel_class'], hover_name=df['ik_2d'])

#rotenone: JUVIOZPCNVVQFO
#deguelin: ORDAZKGHSNRHTD

categories = ['count_in_group', 'specificity', 'sel_pathway', 'sel_class']

dic_cat= {}
for cat in categories:
    if cat == 'count_in_group':
        categorical = False
        colorsIdx = 'inferno_r'
        title = 'Count in group'
    elif cat == 'specificity':
        categorical = False
        colorsIdx = 'inferno_r'
        title = 'Group specificity'
    elif cat == 'sel_pathway':
        categorical = True
        colorsIdx = {
        'npc_Terpenoids': '#ae2012',
        'npc_Shikimates_and_Phenylpropanoids': '#ee9b00',
        'npc_Fatty_acids': '#0a9396',
        'Other': '#e9d8a6'
        }
        title = 'Chemical pathway'
    elif cat == 'sel_class':
      categorical = True
      colorsIdx = {
      'npc_Rotenoids': '#ae2012',
      'npc_Isoflavones': '#ee9b00',
      'Other': '#e9d8a6'
      }
      title = 'Chemical class'
        
    dic_cat[cat] = {}
    dic_cat[cat]['categorical'] = categorical
    dic_cat[cat]['colorsIdx'] = colorsIdx
    dic_cat[cat]['title'] = title
  
  
# create size vector
df['size'] = np.where((df['np_class'] == 'Rotenoids'), 100, 40)  
    
fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, sharey=True, constrained_layout=True)
fig.set_size_inches(30, 9)
fig.set_dpi(300)

for ax, cat in zip([ax1, ax2, ax3, ax4], categories) :
  if cat == 'sel_pathway':
    colors = list(map(dic_cat[cat]['colorsIdx'].get, list(df[cat])))
  else:
    colors = list(df[cat])
  for i in range(len(s)):
    ax.plot(
        [x[s[i]], x[t[i]]],
        [y[s[i]], y[t[i]]],
        "k-",
        alpha=0.5,
        zorder=1,    
    )
    
  if (cat == 'sel_pathway') | (cat == 'sel_class') :
    for pathway in dic_cat[cat]['colorsIdx'].keys():
      df_sel = df[df[cat] == pathway]
      color = dic_cat[cat]['colorsIdx'][pathway]
      ax.set_title(dic_cat[cat]['title'], fontsize=30)     
      scatter = ax.scatter(df_sel['x'], df_sel['y'], zorder=2, s=df_sel['size'], c=color, label=pathway,
                           alpha = 0.8, edgecolors='grey')
      if cat =='sel_pathway':
        legend = ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.3),
            ncol=1, fancybox=False, fontsize=22)
      elif cat =='sel_class':
        legend = ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.25),
            ncol=1, fancybox=False, fontsize=22)
      
  else:
    cmap = dic_cat[cat]['colorsIdx']
    ax.set_title(dic_cat[cat]['title'], fontsize=30)      
    scatter = ax.scatter(x, y, zorder=2, c=colors, cmap=cmap, alpha = 0.8,
                         edgecolors='grey', s=df['size'])
    cbar  = fig.colorbar(scatter, ax=ax, location='bottom', shrink=0.6)
    cbar.ax.tick_params(labelsize=15)
    cbar.ax.set_title(dic_cat[cat]['title'], fontsize=22)
  ax.axis('off')
  
plt.savefig('../output/rotenoids/active_cluster_annotations_tmap.png', bbox_inches='tight')


### CANOPUS Analysis

sparql.setQuery("""
    PREFIX enpkg: <https://enpkg.commons-lab.org/kg/>
    PREFIX enpkgmodule: <https://enpkg.commons-lab.org/module/>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    
    SELECT ?sample ?tcact ?tox ?ccount
    WHERE
    {
      ?sample rdf:type enpkg:LabExtract .
          ?sample enpkgmodule:has_bioassay_results ?tcres .
          ?sample enpkgmodule:has_bioassay_results ?tcres2 .
            ?tcres rdf:type enpkgmodule:Tcruzi10ugml .
            ?tcres2 rdf:type enpkgmodule:L610ugml .
              ?tcres enpkgmodule:inhibition_percentage ?tcact .
              ?tcres2 enpkgmodule:inhibition_percentage ?tox .
      {
    SELECT ?sample (COUNT(DISTINCT ?canopus) AS ?ccount)
    WHERE
    {   
        ?canopus rdf:type enpkg:SiriusCanopusAnnotation .
          ?canopus enpkg:has_canopus_npc_class ?np_class .
          ?canopus enpkg:has_canopus_npc_class_prob ?prob .
        FILTER((regex(str(?np_class), "npc_Rotenoids")) && (?prob > 0.5))
        ?feature enpkg:has_canopus_annotation ?canopus .
          ?feature_list enpkg:has_lcms_feature ?feature .
            ?lcms enpkg:has_lcms_feature_list ?feature_list .
                ?sample enpkg:has_LCMS ?lcms .
    }
    GROUP BY ?sample ORDER BY DESC(?ccount)
       }
    }
""")

results = sparql.queryAndConvert()
df = json_normalize(results['results']["bindings"])
df = df.stack().str.replace(uri_ref, "").unstack()
df.drop(list(df.filter(regex = 'type')), axis = 1, inplace = True)
df.columns = df.columns.str.replace('.value', '')
df['tcact'] = df['tcact'].astype(float)
df['tcact'] = np.where((df['tcact'] >= 100), 100, df['tcact'])
df['tox'] = df['tox'].astype(float)
df['tox'] = np.where((df['tox'] >= 100), 100, df['tox'])

df['ccount'] = df['ccount'].astype(int)
df['hit_status'] = np.where(((df['tcact'] >= 80) & (df['tox'] <= 50)), '#06d6a0', '#e9d8a6')

fig, ax = plt.subplots()
fig.set_size_inches(8, 5)
fig.set_dpi(300)
plt.scatter(x= df['tcact'], y = df['ccount'], c=df['hit_status'], edgecolors='grey')
plt.xlim([-1, df['tcact'].max()+1])
plt.ylim([0, df['ccount'].max()+1])
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
# Only show ticks on the left and bottom spines
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
ax.tick_params(axis='x', labelsize=10)
ax.tick_params(axis='y', labelsize=10)
ax.set_xlabel('Anti-$\it{T. cruzi}$ activity', fontsize=13)
ax.set_ylabel('Number of features annotated\nas rotenoids (CANOPUS, p > 0.5)', fontsize=13)
plt.savefig('../output/rotenoids/canopus_count_rotenoids_vs_activity.png', bbox_inches='tight')



# CANOPUS Analysis ALL CLASSES

selection = 'VGF150_H08|VGF146_F09|VGF146_E08|VGF150_A09|VGF155_D05|VGF139_H07'
#selection = 'VGF141_B11'

sparql.setQuery("""
    PREFIX enpkg: <https://enpkg.commons-lab.org/kg/>
    PREFIX enpkgmodule: <https://enpkg.commons-lab.org/module/>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    
    SELECT * {
      {
        SELECT ?np_class (AVG(?count_sel) AS ?avg_sel)
        WHERE
        {   
          SELECT ?np_class (COUNT(DISTINCT ?canopus) AS ?count_sel)
          WHERE
          {
            ?canopus rdf:type enpkg:SiriusCanopusAnnotation .
              ?canopus enpkg:has_canopus_npc_class ?np_class .
              ?canopus enpkg:has_canopus_npc_class_prob ?prob .
            FILTER((?prob > 0.5))
            ?feature enpkg:has_canopus_annotation ?canopus .
              ?feature_list enpkg:has_lcms_feature ?feature .
                ?lcms enpkg:has_lcms_feature_list ?feature_list .
                  ?sample enpkg:has_LCMS ?lcms .
                FILTER(regex(str(?sample), "%s"))
          } GROUP BY ?sample ?np_class
        } GROUP BY ?np_class
      }
      {
        SELECT ?np_class (AVG(?count_tot) AS ?avg_tot)
        WHERE
        {   
          SELECT ?np_class (COUNT(DISTINCT ?canopus) AS ?count_tot)
          WHERE
          {   
            ?canopus rdf:type enpkg:SiriusCanopusAnnotation .
              ?canopus enpkg:has_canopus_npc_class ?np_class .
              ?canopus enpkg:has_canopus_npc_class_prob ?prob .
            FILTER((?prob > 0.5))
            ?feature enpkg:has_canopus_annotation ?canopus .
              ?feature_list enpkg:has_lcms_feature ?feature .
                ?lcms enpkg:has_lcms_feature_list ?feature_list .
                  ?sample enpkg:has_LCMS ?lcms .
          } GROUP BY ?sample ?np_class
        } GROUP BY ?np_class
      }
    }
""" % (selection)
) 

results = sparql.queryAndConvert()
df2 = json_normalize(results['results']["bindings"])
df2 = df2.stack().str.replace(uri_ref, "").unstack()
df2.drop(list(df2.filter(regex = 'type')), axis = 1, inplace = True)
df2.columns = df2.columns.str.replace('.value', '')
df2['avg_sel'] = df2['avg_sel'].astype(float)
df2['avg_tot'] = df2['avg_tot'].astype(float)

df2['selected_class'] = np.where(df2['np_class'].isin(['npc_Rotenoids', 'npc_Isoflavones', 'npc_Flavanones', 'npc_Other_Octadecanoids', 'npc_Diacylglycerols']), \
  df2['np_class'], 'Other')

color_map = {
  'npc_Rotenoids': '#ae2012',
  'npc_Isoflavones':'#ee9b00',
  'npc_Flavanones':'#ca6702',
  'npc_Other_Octadecanoids': '#94d2bd',
  'npc_Diacylglycerols': '#001219',
  'Other': '#e9d8a6'
}

colors = df2['selected_class'].map(color_map)
    
fig, ax = plt.subplots()
fig.set_size_inches(8, 5)
fig.set_dpi(600)

for np_class in df2['selected_class'].unique():
  submerge = df2[df2['selected_class'] == np_class].copy()
  plt.scatter(x= submerge['avg_tot'], y = submerge['avg_sel'], c=color_map[np_class], edgecolors='grey', label = np_class)

legend = ax.legend(loc="upper right", ncol=1, fancybox=False, fontsize=11)
# ax.plot([0,1],[0,1], c = 'grey')

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
# Only show ticks on the left and bottom spines
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
ax.tick_params(axis='x', labelsize=10)
ax.tick_params(axis='y', labelsize=10)
ax.set_ylabel('Average CANOPUS annotation count (p > 0.5)\nin active selection (n=6)', fontsize=11)
ax.set_xlabel('Average CANOPUS annotation count (p > 0.5)\nin the whole set (n=1,600)', fontsize=11)
xpoints = ypoints = ax.get_xlim()
ax.plot(xpoints, ypoints, color='grey', scalex=False, scaley=False)
plt.xlim([0, df2['avg_tot'].max()+1])
plt.ylim([0, df2['avg_sel'].max()+1])

angle = 45
ax.text(*np.array((12.6, 9.5)), 'Under-annotated', fontsize=13, c='grey',
          rotation=angle, rotation_mode='anchor',
          transform_rotates_text=True)
ax.text(*np.array((12.5, 14)), 'Over-annotated', fontsize=13, c='grey',
          rotation=angle, rotation_mode='anchor',
          transform_rotates_text=True)

plt.savefig('../output/rotenoids/canopus_all_classes_difference.png', bbox_inches='tight', dpi=300)

px.scatter(x= df2['avg_tot'], y = df2['avg_sel'], color=df2['np_class'])