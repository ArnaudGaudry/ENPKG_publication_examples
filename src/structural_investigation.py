import pandas as pd
from pandas import json_normalize
from SPARQLWrapper import SPARQLWrapper, JSON
from rdkit.Chem import PandasTools
import pandas as pd
import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from tqdm import tqdm
import plotly.express as px
import tmap as tm
from mhfp.encoder import MHFPEncoder
from faerun import Faerun
import numpy as np
import plotly.graph_objects as go


# Fonctions

def compute_tmap(table):
    
    enc = MHFPEncoder(1024)
    lf = tm.LSHForest(1024, 64)
    
    i = 0
    fps = []
    for _, row in table.iterrows():
        if i != 0 and i % 10 == 0:
                print(100 * i / len(table))
        mol = AllChem.MolFromSmiles(row["smiles"])
        fps.append(tm.VectorUint(enc.encode_mol(mol)))        
        i+=1

    lf.batch_add(fps)
    lf.index()

    CFG_TMAP = tm.LayoutConfiguration()
    CFG_TMAP.k = 50
    CFG_TMAP.kc = 50

    x, y, s, t, _ = tm.layout_from_lsh_forest(lf, CFG_TMAP)
    return x, y, s, t

# SPARQL
uri_ref = 'https://enpkg.commons-lab.org/kg/'
uri_module = 'https://enpkg.commons-lab.org/module/'

sparql = SPARQLWrapper('http://FASIE-1439:7200/repositories/graph_sandbox')
sparql.setReturnFormat(JSON)

# Get all ChEMBL annotation with an activity given as an IC50 against L. donovani

sparql.setQuery("""
    PREFIX enpkg: <https://enpkg.commons-lab.org/kg/>
    PREFIX enpkgmodule: <https://enpkg.commons-lab.org/module/>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>

    SELECT ?ik (SAMPLE(?smiles) AS ?smiles) (SAMPLE(?activity) AS ?activity) (SAMPLE(?feature) AS ?feature)
    WHERE
    {   
        ?ik rdf:type enpkg:InChIkey .
        ?ik enpkgmodule:has_chembl_id ?chembl_id .
        ?ik enpkg:has_smiles ?smiles .
        ?chembl_id enpkgmodule:has_chembl_activity ?chembl_activity .
            ?chembl_activity enpkgmodule:target_id ?target .
            ?chembl_activity enpkgmodule:activity_type ?act_type .
            ?chembl_activity enpkgmodule:activity_value ?activity .
        FILTER(regex(str(?target), "CHEMBL367"))
        FILTER(regex(str(?act_type), "IC50"))
        FILTER(?activity < 500)
    OPTIONAL
        {?ik2d enpkg:is_InChIkey2D_of ?ik .
            ?siriusannotation enpkg:has_InChIkey2D ?ik2d .
            ?isdbannotation enpkg:has_InChIkey2D ?ik2d . 
                ?feature enpkg:has_sirius_annotation ?siriusannotation .
                ?feature enpkg:has_isdb_annotation ?isdbannotation .
                {?siriusannotation rdf:type enpkg:SiriusStructureAnnotation .
                    ?siriusannotation enpkg:has_cosmic_score ?cosmic .
                    ?siriusannotation enpkg:has_zodiac_score ?zodiac .
                    FILTER((?cosmic > 0.5) && (?zodiac > 0.8))}
                UNION
                {?isdbannotation rdf:type enpkg:IsdbAnnotation .
                    ?isdbannotation enpkg:has_isdb_taxo_score ?taxo_score .
                    FILTER(?taxo_score >= 6)}
            }
    }
    GROUP BY ?ik
""")

results = sparql.queryAndConvert()
df = json_normalize(results['results']["bindings"])
df = df.stack().str.replace(uri_ref, "").unstack()
df.drop(list(df.filter(regex = 'type')), axis = 1, inplace = True)
df.columns = df.columns.str.replace('.value', '')
PandasTools.AddMoleculeColumnToFrame(df,'smiles','Molecule')

df['feature'] = df['feature'].fillna('none')
df['feature'] = np.where((df['feature'] != 'none'), 'Yes', 'No')
df.sort_values(by=['feature'], axis=0, inplace=True, ascending=True)

  
x, y, s, t = compute_tmap(df)

edges_list = []
for source, target in zip(s,t):
    x_source = x[source]
    y_source = y[source]
    x_target = x[target]
    y_target = y[target]
    coord = ((x_source, y_source), (x_target, y_target))
    edges_list.append(coord)
    
df['x'] = x
df['y'] = y
df['activity']=df['activity'].astype(float)
df['size'] = np.where((df['feature'] == 'Yes'), 10, 5)


categories = ['feature', 'activity']
dic_cat= {}
for cat in categories:
    if cat == 'activity':
        categorical = False
        colorsIdx = 'Inferno'
        title = 'Activity'
    elif cat == 'feature':
        categorical = True
        colorsIdx = {
            'Yes': '#ae2012',
            'No': '#e9d8a6'
        }
        title = 'Annotation'
    dic_cat[cat] = {}
    dic_cat[cat]['categorical'] = categorical
    dic_cat[cat]['colorsIdx'] = colorsIdx
    dic_cat[cat]['title'] = title


# fig = px.scatter(x = df['x'], y = df['y'], color = df['feature'],
#                  color_discrete_sequence=["#ae2012", "#e9d8a6"],
#                  hover_name=df['ik'])
# fig.update_layout(template='simple_white')

category = 'activity'

fig = go.Figure()
# Plot tree
for edge in edges_list:
    from_x = edge[0][0]
    from_y = edge[0][1]
    to_x = edge[1][0]
    to_y = edge[1][1]
    fig.add_trace(
        go.Scatter(
            x= [from_x, to_x],
            y = [from_y, to_y],
            mode='lines', showlegend = False,
            line=dict(color='black', width=0.07)
        ))

cats = list(df[category].unique())

if dic_cat[category]['categorical'] is True:
    for cat in cats:  
        result_cat = df[df[category] == cat]
        if (cat == 'No'):
            size=7
        else:
            size=7
    
        fig.add_trace(
            go.Scatter(
                x= result_cat['x'], y = result_cat['y'],
                mode='markers',
                marker=dict(
                    size = size,
                    opacity=0.85,
                    line_width=0.5,
                    line_color ='grey'
                ),
                marker_color= dic_cat[category]['colorsIdx'][cat], name=cat, legendgroup=cat,
                hovertext = result_cat['ik'], line_width=1,
                line_color ='grey'
                )
            )

elif dic_cat[category]['categorical'] is False:
    fig.add_trace(
        go.Scatter(
            opacity=1,
            mode='markers',
            x= df['x'], y = df['y'],
            marker=dict(
                size = 7,
                line_width=0.5,
                line_color ='lightgray',
                color=df[category], 
                colorscale='Inferno_r', 
                showscale=True
            ),
            name=category,
            hovertext = df['ik']
            ),
        )


fig.update_layout(height=600, width=600, template = 'simple_white')

fig.update_layout(
    legend=dict(
        title = dic_cat[category]['title'],
        orientation="h",
        y=-0,
        x = 0,
        bordercolor="Black",
        borderwidth=2
        ),
    font=dict(
        family="Arial",
        size=16,
        color="black"
        )    
    )

fig.update_annotations(font_size=20)
fig.update_xaxes(visible=False)
fig.update_yaxes(visible=False)

fig.show()
  
fig.write_html(f"/mnt/c/Users/gaudrya.FARMA/Desktop/chembl_annot_color_{category}.html")
fig.write_image(f"/mnt/c/Users/gaudrya.FARMA/Desktop/chembl_annot_color_{category}.jpeg",  scale=3)
fig.write_image(f"/mnt/c/Users/gaudrya.FARMA/Desktop/chembl_annot_color_{category}.svg",  scale=3)



# Sachem similarity

# Query to return Sirius annotations in a given samples, ordered by relative intensity
sparql.setQuery("""
    PREFIX enpkg: <https://enpkg.commons-lab.org/kg/>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    PREFIX wdt: <http://www.wikidata.org/prop/direct/>
    PREFIX sachem: <http://bioinfo.uochb.cas.cz/rdf/v1.0/sachem#> # prefixes needed for structural similarity search
    PREFIX idsm: <https://idsm.elixir-czech.cz/sparql/endpoint/>
    PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
    
    SELECT ?ik2d (SAMPLE(?smiles) AS ?smiles) (SAMPLE(?SCORE) AS ?SCORE) #?feature ?annotation ?ik2d ?rt ?rel_area ?smiles
    WHERE
    {   ?sample rdf:type enpkg:LabExtract .
        FILTER(regex(str(?sample), "VGF141_B11"))
          ?sample enpkg:has_LCMS ?lcms .
            ?lcms rdf:type enpkg:LCMSAnalysisPos .
            ?lcms enpkg:has_lcms_feature_list ?feature_list .
              ?feature_list enpkg:has_lcms_feature ?feature .
                {?feature enpkg:has_sirius_annotation ?annotation .
                  ?annotation enpkg:has_InChIkey2D ?ik2d .
                  ?annotation enpkg:has_cosmic_score ?cosmic .
                  ?annotation enpkg:has_zodiac_score ?zodiac .
                  FILTER((?cosmic > 0.5) && (?zodiac > 0.8))}
                  UNION
                {?feature enpkg:has_isdb_annotation ?annotation .
                  ?annotation enpkg:has_taxo_score ?taxo .
                  FILTER(?taxo >= 6)}
                  
                  ?annotation enpkg:has_InChIkey2D ?ik2d .
                    ?ik2d enpkg:has_smiles ?smiles .
                    ?ik2d enpkg:is_InChIkey2D_of ?ik .
                      ?ik enpkg:has_wd_id ?wd_id .
                      SERVICE idsm:wikidata {
                        SELECT * WHERE {
                            [ sachem:compound ?wd_id;
                                sachem:score ?SCORE ] sachem:similaritySearch [
                                    sachem:query "CC1=C(C(=O)C=C2C1=CC=C3C2(CCC4(C3(CCC5(C4CC(=C)CC5)C)C)C)C)O";
                                    sachem:cutoff "0.01"^^xsd:double ].
                      }   }                  
                ?feature enpkg:has_retention_time ?rt .
                ?feature enpkg:has_relative_feature_area ?rel_area .
    } GROUP BY ?ik2d
""")


results = sparql.queryAndConvert()
df = json_normalize(results['results']["bindings"])
df = df.stack().str.replace(uri_ref, "").unstack()
df.drop(list(df.filter(regex = 'type')), axis = 1, inplace = True)
df.columns = df.columns.str.replace('.value', '')
df.drop_duplicates(inplace=True)
PandasTools.AddMoleculeColumnToFrame(df,'smiles','Molecule')

x, y, s, t = compute_tmap(df)

px.scatter(x= x, y = y)

edges_list = []
for source, target in zip(s,t):
    x_source = x[source]
    y_source = y[source]
    x_target = x[target]
    y_target = y[target]
    coord = ((x_source, y_source), (x_target, y_target))
    edges_list.append(coord)
    
df['x'] = x
df['y'] = y
df['SCORE']=df['SCORE'].astype(float)
df['cutoff'] = np.where((df['SCORE'] > 0.8), 'gold', 'gray')

fig = go.Figure()
# Plot tree
for edge in edges_list:
    from_x = edge[0][0]
    from_y = edge[0][1]
    to_x = edge[1][0]
    to_y = edge[1][1]
    fig.add_trace(
        go.Scatter(
            x= [from_x, to_x],
            y = [from_y, to_y],
            mode='lines', showlegend = False,
            line=dict(color='black', width=0.07)
        ))

fig.add_trace(
    go.Scatter(
        opacity=1,
        mode='markers',
        x= df['x'], y = df['y'],
        marker=dict(
            size = 10,
            line_width=2,
            line_color =df['cutoff'],
            color=df['SCORE'], 
            colorscale='Inferno_r', 
            showscale=True
        ),
        #name = category,
        hovertext = df['ik2d']
        ),
    )


fig.update_layout(height=600, width=600, template = 'simple_white')

fig.update_layout(
    legend=dict(
        #title = dic_cat[category]['title'],
        orientation="h",
        y=-0,
        x = 0,
        bordercolor="Black",
        borderwidth=2
        ),
    font=dict(
        family="Arial",
        size=16,
        color="black"
        )    
    )

fig.update_annotations(font_size=20)
fig.update_xaxes(visible=False)
fig.update_yaxes(visible=False)

fig.show()

  
fig.write_html(f"/mnt/c/Users/gaudrya.FARMA/Desktop/pristimera_color_similarity.html")
fig.write_image(f"/mnt/c/Users/gaudrya.FARMA/Desktop/pristimera_color_similarity.jpeg",  scale=3)
fig.write_image(f"/mnt/c/Users/gaudrya.FARMA/Desktop/pristimera_color_similarity.svg",  scale=3)



  