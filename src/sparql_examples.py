import pandas as pd
import numpy as np
from pandas import json_normalize
from SPARQLWrapper import SPARQLWrapper, JSON
from rdkit.Chem import PandasTools
import plotly.express as px

uri_ref = 'https://enpkg.commons-lab.org/kg/'
uri_module = 'https://enpkg.commons-lab.org/module/'

sparql = SPARQLWrapper('https://enpkg.commons-lab.org/graphdb/repositories/ENPKG')
#sparql = SPARQLWrapper('http://FASIE-1439:7200/repositories/graph_sandbox')
sparql.setReturnFormat(JSON)


# Query to return the number of features in PI mode

sparql.setQuery("""
    PREFIX enpkg: <https://enpkg.commons-lab.org/kg/>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>

    SELECT (COUNT(?feature) as ?count)
    WHERE
    { 
        ?feature rdf:type enpkg:LCMSFeature .
        ?featurelist enpkg:has_lcms_feature ?feature .
        ?lcms enpkg:has_lcms_feature_list ?featurelist .
        ?lcms rdf:type enpkg:LCMSAnalysisPos
    }
""")

results = sparql.queryAndConvert()

df = json_normalize(results['results']["bindings"])
df = df.stack().str.replace(uri_ref, "").str.replace(uri_module, "").unstack()
df.drop(list(df.filter(regex = 'type')), axis = 1, inplace = True)
df.columns = df.columns.str.replace('.value', '')
df


# Query to return the number of features in NI mode

sparql.setQuery("""
    PREFIX enpkg: <https://enpkg.commons-lab.org/kg/>
    PREFIX enpkgmodule: <https://enpkg.commons-lab.org/module/>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    PREFIX wdt: <http://www.wikidata.org/prop/direct/>
    PREFIX wd: <http://www.wikidata.org/entity/>

    SELECT (COUNT(?feature) as ?count)
    WHERE
    { 
        ?feature rdf:type enpkg:LCMSFeature .
        ?featurelist enpkg:has_lcms_feature ?feature .
        ?lcms enpkg:has_lcms_feature_list ?featurelist .
        ?lcms rdf:type enpkg:LCMSAnalysisNeg
    }
""")

results = sparql.queryAndConvert()

df = json_normalize(results['results']["bindings"])
df = df.stack().str.replace(uri_ref, "").str.replace(uri_module, "").unstack()
df.drop(list(df.filter(regex = 'type')), axis = 1, inplace = True)
df.columns = df.columns.str.replace('.value', '')
df

# Query to return the number  of features with CANOPUS annotations

sparql.setQuery("""
    PREFIX enpkg: <https://enpkg.commons-lab.org/kg/>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>

    SELECT (COUNT(DISTINCT(?feature)) as ?count)
    WHERE
    {   ?feature rdf:type enpkg:LCMSFeature .
        ?feature enpkg:has_canopus_annotation ?annotation .
    }
""")

results = sparql.queryAndConvert()

df = json_normalize(results['results']["bindings"])
df = df.stack().str.replace(uri_ref, "").str.replace(uri_module, "").unstack()
df.drop(list(df.filter(regex = 'type')), axis = 1, inplace = True)
df.columns = df.columns.str.replace('.value', '')
df

# Query to return the number of features with SIRIUS annotations

sparql.setQuery("""
    PREFIX enpkg: <https://enpkg.commons-lab.org/kg/>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>

    SELECT (COUNT(DISTINCT(?feature)) as ?count)
    WHERE
    {   ?feature rdf:type enpkg:LCMSFeature .
        ?feature enpkg:has_sirius_annotation ?annotation .
    }
""")
results = sparql.queryAndConvert()

df = json_normalize(results['results']["bindings"])
df = df.stack().str.replace(uri_ref, "").str.replace(uri_module, "").unstack()
df.drop(list(df.filter(regex = 'type')), axis = 1, inplace = True)
df.columns = df.columns.str.replace('.value', '')
df

# Query to return the number  of features with ISDB annotations

sparql.setQuery("""
    PREFIX enpkg: <https://enpkg.commons-lab.org/kg/>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>

    SELECT (COUNT(DISTINCT(?feature)) as ?count)
    WHERE
    {   ?feature rdf:type enpkg:LCMSFeature .
        ?feature enpkg:has_isdb_annotation ?annotation .
    }
""")
results = sparql.queryAndConvert()

df = json_normalize(results['results']["bindings"])
df = df.stack().str.replace(uri_ref, "").str.replace(uri_module, "").unstack()
df.drop(list(df.filter(regex = 'type')), axis = 1, inplace = True)
df.columns = df.columns.str.replace('.value', '')
df

# Query to return the number distinct IK2D in annotations

sparql.setQuery("""
    PREFIX enpkg: <https://enpkg.commons-lab.org/kg/>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>

    SELECT (COUNT(DISTINCT(?ik2d)) as ?count)
    WHERE
    {   ?annotation rdf:type enpkg:Annotation .
        ?annotation enpkg:has_InChIkey2D ?ik2d .
    }
""")
results = sparql.queryAndConvert()

df = json_normalize(results['results']["bindings"])
df = df.stack().str.replace(uri_ref, "").str.replace(uri_module, "").unstack()
df.drop(list(df.filter(regex = 'type')), axis = 1, inplace = True)
df.columns = df.columns.str.replace('.value', '')
df

# Query to return the number of wd ID from IK2D in annotations

sparql.setQuery("""
    PREFIX enpkg: <https://enpkg.commons-lab.org/kg/>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>

    SELECT (COUNT(DISTINCT(?wd)) as ?count)
    WHERE
    {   ?annotation rdf:type enpkg:Annotation .
        ?annotation enpkg:has_InChIkey2D ?ik2d .
        ?ik2d enpkg:is_InChIkey2D_of ?ik .
        ?ik enpkg:has_wd_id ?wd
    }
""")
results = sparql.queryAndConvert()

df = json_normalize(results['results']["bindings"])
df = df.stack().str.replace(uri_ref, "").str.replace(uri_module, "").unstack()
df.drop(list(df.filter(regex = 'type')), axis = 1, inplace = True)
df.columns = df.columns.str.replace('.value', '')
df

# Query to return samples active against T cruzi without cytotoxicity

sparql.setQuery("""
    PREFIX enpkg: <https://enpkg.commons-lab.org/kg/>
    PREFIX enpkgmodule: <https://enpkg.commons-lab.org/module/>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    PREFIX wdt: <http://www.wikidata.org/prop/direct/>
    PREFIX wd: <http://www.wikidata.org/entity/>

    SELECT ?active_extract ?material ?wd_id ?taxon ?organe
    WHERE
        { 
            service <https://query.wikidata.org/sparql> {
                    ?wd_id wdt:P225 ?taxon .
                }
            { SELECT ?active_extract ?material ?wd_id ?organe WHERE
        { 
            ?active_extract rdf:type enpkg:LabExtract .
                ?active_extract enpkgmodule:has_bioassay_results ?biores .
                ?active_extract enpkgmodule:has_bioassay_results ?toxres .
                    ?biores rdf:type enpkgmodule:Tcruzi10ugml .
                    ?toxres rdf:type enpkgmodule:L610ugml .
                        ?biores enpkgmodule:inhibition_percentage ?tc_inhib .
                        ?toxres enpkgmodule:inhibition_percentage ?l6_inhib .
                        FILTER((?tc_inhib > 80) && (?l6_inhib < 50))
            ?material enpkg:has_lab_process ?active_extract .
            ?material enpkgmodule:has_organe ?organe .
                ?material enpkg:has_wd_id ?wd_id .
                }
            }
        }
""")

results = sparql.queryAndConvert()

df = json_normalize(results['results']["bindings"])
df = df.stack().str.replace(uri_ref, "").str.replace(uri_module, "").unstack()
df.drop(list(df.filter(regex = 'type')), axis = 1, inplace = True)
df.columns = df.columns.str.replace('.value', '')
df


# Query 1: Return the number of features with the same sirius and isdb annotation

sparql.setQuery("""
    PREFIX enpkg: <https://enpkg.commons-lab.org/kg/>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>

    SELECT (COUNT(?feature) AS ?count) 
    WHERE
    {   
        ?feature rdf:type enpkg:LCMSFeature .
          ?feature enpkg:has_sirius_annotation ?sirius_annotation .
            ?sirius_annotation rdf:type enpkg:SiriusStructureAnnotation .
              ?sirius_annotation enpkg:has_InChIkey2D ?sirius_ik2d .
          ?feature enpkg:has_isdb_annotation ?isdb_annotation .
            ?isdb_annotation rdf:type enpkg:IsdbAnnotation .
              ?isdb_annotation enpkg:has_InChIkey2D ?isdb_ik2d .
        FILTER(?isdb_ik2d = ?sirius_ik2d) .
    }
""")

results = sparql.queryAndConvert()
df = json_normalize(results['results']["bindings"])
df = df.stack().str.replace(uri_ref, "").str.replace(uri_module, "").unstack()
df.drop(list(df.filter(regex = 'type')), axis = 1, inplace = True)
df.columns = df.columns.str.replace('.value', '')
df


# Query 2: Return the samples withe features annotated as Aspidosperma type alkaloids by Canopus

sparql.setQuery("""
    PREFIX enpkg: <https://enpkg.commons-lab.org/kg/>
    PREFIX enpkgmodule: <https://enpkg.commons-lab.org/module/>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    PREFIX wdt: <http://www.wikidata.org/prop/direct/>
    PREFIX wd: <http://www.wikidata.org/entity/>

    SELECT ?extract ?organe ?species_name ?genus_name ?family_name ?count_of_selected_class
    WHERE
     {  
        ?material enpkg:has_lab_process ?extract .
            ?material enpkgmodule:has_organe ?organe .
            ?material enpkg:has_wd_id ?wd_sp .
            OPTIONAL{
                SERVICE <https://query.wikidata.org/sparql> {
                ?wd_sp wdt:P225 ?species_name .
                ?family wdt:P31 wd:Q16521 ;
                    wdt:P105 wd:Q35409 ;
                    wdt:P225 ?family_name ;
                    ^wdt:P171* ?wd_sp .
                ?genus wdt:P31 wd:Q16521 ;
                    wdt:P105 wd:Q34740 ;
                    wdt:P225 ?genus_name ;
                    ^wdt:P171* ?wd_sp 
                }
            }
        {
            SELECT ?extract (COUNT(DISTINCT ?feature) AS ?count_of_selected_class)
            WHERE
            {   
                ?extract rdf:type enpkg:LabExtract .
                    ?extract enpkg:has_LCMS ?lcms .
                        ?lcms rdf:type ?LCMSAnalysisPos .
                        ?lcms enpkg:has_lcms_feature_list ?feature_list .
                            ?feature_list enpkg:has_lcms_feature ?feature .
                                ?feature enpkg:has_canopus_annotation ?canopus .
                                    ?canopus enpkg:has_canopus_npc_class ?np_class .
                                    ?canopus enpkg:has_canopus_npc_class_prob ?class_prob .
                                    FILTER(regex(str(?np_class), "Aspidosperma_type")) .
                                    FILTER((?class_prob > 0.5)) .
            } GROUP BY ?extract ORDER BY DESC(?count_of_selected_class)
        }
    }

""")

results = sparql.queryAndConvert()
df = json_normalize(results['results']["bindings"])
df = df.stack().str.replace(uri_ref, "").unstack()
df.drop(list(df.filter(regex = 'type')), axis = 1, inplace = True)
df.columns = df.columns.str.replace('.value', '')
df

# Query 3: Get features where annotations contain a given substructure

sparql.setQuery("""
    PREFIX enpkg: <https://enpkg.commons-lab.org/kg/>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    PREFIX wdt: <http://www.wikidata.org/prop/direct/>
    PREFIX wd: <http://www.wikidata.org/entity/>
    PREFIX idsm: <https://idsm.elixir-czech.cz/sparql/endpoint/>
    PREFIX sachem: <http://bioinfo.uochb.cas.cz/rdf/v1.0/sachem#>

    SELECT DISTINCT ?ik2d ?smiles
    WHERE { 
        ?extract rdf:type enpkg:LabExtract .
        FILTER(regex(str(?extract), "VGF152_B02"))
            ?extract enpkg:has_LCMS ?lcms .
                ?lcms enpkg:has_lcms_feature_list ?feature_list .
                ?feature_list enpkg:has_lcms_feature ?feature .
                    ?feature enpkg:has_sirius_annotation|enpkg:has_isdb_annotation ?annotation . 
                    ?annotation enpkg:has_InChIkey2D ?ik2d .
                        ?ik2d enpkg:has_smiles ?smiles .
                        ?ik2d enpkg:is_InChIkey2D_of ?ik .
                            ?ik enpkg:has_wd_id ?wd_id .
                            SERVICE idsm:wikidata {
                                VALUES ?SUBSTRUCTURE {
                                "CCC12CCCN3C1C4(CC3)C(CC2)NC5=CC=CC=C45" # Aspidospermidine scaffold
                                }
                                ?wd_id sachem:substructureSearch _:b16.
                                _:b16 sachem:query ?SUBSTRUCTURE.
                            }      
    }
""")

results = sparql.queryAndConvert()
df = json_normalize(results['results']["bindings"])
df = df.stack().str.replace(uri_ref, "").unstack()
df.drop(list(df.filter(regex = 'type')), axis = 1, inplace = True)
df.columns = df.columns.str.replace('.value', '')
df


# Query 5: Get structural annotations of an extract reported in the same genus

sparql.setQuery("""
    PREFIX enpkg: <https://enpkg.commons-lab.org/kg/>
    PREFIX enpkgmodule: <https://enpkg.commons-lab.org/module/>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    PREFIX wdt: <http://www.wikidata.org/prop/direct/>
    PREFIX wd: <http://www.wikidata.org/entity/>
    PREFIX prov: <http://www.w3.org/ns/prov#>
    PREFIX pr: <http://www.wikidata.org/prop/reference/>
    
    SELECT DISTINCT ?ik2d ?genus WHERE
    {   ?material enpkg:has_lab_process ?extract .
            ?material   enpkgmodule:has_organe ?organe ;
                        enpkg:has_wd_id ?wd_sp .
                FILTER(regex(str(?extract), "VGF152_B02"))
                ?extract enpkg:has_LCMS ?lcms .
                    ?lcms rdf:type ?LCMSAnalysisPos .
                    ?lcms enpkg:has_lcms_feature_list ?feature_list .
                        ?feature_list enpkg:has_lcms_feature ?feature .
                            ?feature enpkg:has_sirius_annotation ?annotation . 
                                ?annotation enpkg:has_InChIkey2D ?ik2d .
                                    ?ik2d enpkg:has_smiles ?smiles .
                                    ?ik2d enpkg:is_InChIkey2D_of ?ik .
                                        ?ik enpkg:has_wd_id ?wd_id .
        {
        SELECT DISTINCT ?wd_id ?genus WHERE {
            ?material enpkg:has_lab_process ?extract ;
                enpkg:has_wd_id ?wd_sp .
                FILTER(regex(str(?extract), "VGF152_B02")) .
                OPTIONAL{
                    service <https://query.wikidata.org/sparql> {
                        ?wd_sp wdt:P225 ?species_name .
                        ?genus wdt:P31 wd:Q16521 ;
                                wdt:P105 wd:Q34740 ;
                                ^wdt:P171* ?wd_sp .
                        ?childtaxa wdt:P171* ?genus .
                        ?wd_id wdt:P703 ?childtaxa
                    }
                }
            }
    }
}

""")

results = sparql.queryAndConvert()
df = json_normalize(results['results']["bindings"])
df = df.stack().str.replace(uri_ref, "").unstack()
df.drop(list(df.filter(regex = 'type')), axis = 1, inplace = True)
df.columns = df.columns.str.replace('.value', '')
df


# Query 6:  Return compounds annotated in Melochia and active in ChEMBL, \
    # with the taxa they are reported in in Wikidata

sparql.setQuery("""
    PREFIX enpkg: <https://enpkg.commons-lab.org/kg/>
    PREFIX enpkgmodule: <https://enpkg.commons-lab.org/module/>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    PREFIX wdt: <http://www.wikidata.org/prop/direct/>
    PREFIX wd: <http://www.wikidata.org/entity/>

    SELECT DISTINCT ?ik ?chemblid ?name ?type ?value ?wd_id ?taxon ?taxon_name
    WHERE
    { service <https://query.wikidata.org/sparql> {
            ?wd_id wdt:P31 wd:Q11173 .
                ?wd_id wdt:P703 ?taxon .
                ?taxon wdt:P225 ?taxon_name .
        }   
        { SELECT DISTINCT ?ik ?chemblid ?name ?type ?value ?wd_id
            WHERE
            { 
                ?sample rdf:type enpkg:LabExtract
                FILTER(regex(str(?sample), "VGF156_A06"))
                ?sample enpkg:has_LCMS ?lcms .
                    ?lcms rdf:type ?LCMSAnalysisPos .
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
                        ?ik enpkgmodule:has_chembl_id ?chemblid .
                        ?ik enpkg:has_wd_id ?wd_id .
                            ?chemblid enpkgmodule:has_chembl_activity ?chembl_activity .
                            ?chembl_activity enpkgmodule:target_name ?name .
                            FILTER(regex(str(?name), "Trypanosoma cruzi")) 
                            ?chembl_activity enpkgmodule:activity_type ?type .
                            ?chembl_activity enpkgmodule:activity_value ?value .
        }
        }      
    }
""")

results = sparql.queryAndConvert()
df = json_normalize(results['results']["bindings"])
df = df.stack().str.replace(uri_ref, "").unstack()
df.drop(list(df.filter(regex = 'type')), axis = 1, inplace = True)
df.columns = df.columns.str.replace('.value', '')
df


# Query 7:  Query to return features with the highest spec2vec overlap with a given feature

sparql.setQuery("""
    PREFIX enpkg: <https://enpkg.commons-lab.org/kg/>
    PREFIX enpkgmodule: <https://enpkg.commons-lab.org/module/>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    PREFIX wdt: <http://www.wikidata.org/prop/direct/>
              
    SELECT ?feature (SAMPLE(?rt) AS ?rt) (SAMPLE(?parent_mass) AS ?parent_mass) (COUNT(?peakloss) AS ?count) WHERE
          { 
            ?feature rdf:type enpkg:LCMSFeature .
            ?feature enpkg:has_spec2vec_doc ?doc .
            ?feature enpkg:has_parent_mass ?parent_mass .
            ?feature enpkg:has_retention_time ?rt .
              ?doc enpkg:has_spec2vec_loss|enpkg:has_spec2vec_peak ?peakloss .
        
         {SELECT ?peakloss WHERE {
          ?feature rdf:type enpkg:LCMSFeature
          FILTER(regex(str(?feature), "SC_AP_Wi_DCM_features_ms2_pos.mgf:scan:1$"))
            ?feature enpkg:has_spec2vec_doc ?doc .
              ?doc enpkg:has_spec2vec_loss|enpkg:has_spec2vec_peak ?peakloss .
            }
           }
         } GROUP BY ?feature ORDER BY DESC(?count)
""")

results = sparql.queryAndConvert()
df = json_normalize(results['results']["bindings"])
df = df.stack().str.replace(uri_ref, "").unstack()
df.drop(list(df.filter(regex = 'type')), axis = 1, inplace = True)
df.columns = df.columns.str.replace('.value', '')
df

# Query 8: get the features, annotated as M+H by SIRIUS, for which the M-H adduct is detected in the same RT range for a given sample

sparql.setQuery("""
    PREFIX enpkg: <https://enpkg.commons-lab.org/kg/>
    PREFIX enpkgmodule: <https://enpkg.commons-lab.org/module/>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    PREFIX wdt: <http://www.wikidata.org/prop/direct/>
              
    SELECT DISTINCT ?feature ?rt ?pm ?feature_opp ?rt_opp ?pm_opp
    WHERE
    { 
        VALUES ?ppm {
            "5"^^xsd:decimal
            }
        
        ?sample rdf:type enpkg:LabExtract
        FILTER(regex(str(?sample), "VGF156_A06"))
        
        ?sample enpkg:has_LCMS ?lcms .
            ?lcms rdf:type enpkg:LCMSAnalysisPos .
            ?lcms enpkg:has_lcms_feature_list ?feature_list .
                ?feature_list enpkg:has_lcms_feature ?feature .                    
                    ?feature enpkg:has_parent_mass ?pm .
                    ?feature enpkg:has_retention_time ?rt .
                    ?feature enpkg:has_sirius_annotation ?sirius .
                        ?sirius enpkg:has_sirius_adduct ?adduct .
                        FILTER(regex(str(?adduct), "[M+H]+"))
            
        ?sample enpkg:has_LCMS ?lcms_opp .
        ?lcms_opp rdf:type enpkg:LCMSAnalysisNeg .
        ?lcms_opp enpkg:has_lcms_feature_list ?feature_list_opp .
            ?feature_list_opp enpkg:has_lcms_feature ?feature_opp .
            ?feature_opp enpkg:has_parent_mass ?pm_opp .
            ?feature_opp enpkg:has_retention_time ?rt_opp .
        FILTER(((?rt - 0.05) < ?rt_opp) && ((?rt + 0.05) > ?rt_opp))
        FILTER((?pm_opp > ((?pm - 2.014) - ((?ppm * 0.000001) * (?pm - 2.014)))) && (?pm_opp < ((?pm - 2.014) + ((?ppm * 0.000001) * (?pm - 2.014)))))
    }
""")


results = sparql.queryAndConvert()
df = json_normalize(results['results']["bindings"])
df = df.stack().str.replace(uri_ref, "").unstack()
df.drop(list(df.filter(regex = 'type')), axis = 1, inplace = True)
df.columns = df.columns.str.replace('.value', '')
df

# Query 9: get the features wor which a feature with the same annotation is found at the same RT

sparql.setQuery("""
    PREFIX enpkg: <https://enpkg.commons-lab.org/kg/>
    PREFIX enpkgmodule: <https://enpkg.commons-lab.org/module/>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    PREFIX wdt: <http://www.wikidata.org/prop/direct/>
              
    SELECT DISTINCT ?feature ?feature_opp ?ik2d ?rt ?rt_opp
    WHERE
    { 
        ?sample rdf:type enpkg:LabExtract
        FILTER(regex(str(?sample), "VGF156_A06"))
        
        ?sample enpkg:has_LCMS ?lcms .
            ?lcms rdf:type enpkg:LCMSAnalysisPos .
            ?lcms enpkg:has_lcms_feature_list ?feature_list .
                ?feature_list enpkg:has_lcms_feature ?feature .                    
                    ?feature enpkg:has_retention_time ?rt .
                    ?feature enpkg:has_sirius_annotation ?sirius .
                        ?sirius enpkg:has_InChIkey2D ?ik2d .
                        
            
        ?sample enpkg:has_LCMS ?lcms_opp .
        ?lcms_opp rdf:type enpkg:LCMSAnalysisNeg .
        ?lcms_opp enpkg:has_lcms_feature_list ?feature_list_opp .
            ?feature_list_opp enpkg:has_lcms_feature ?feature_opp .
                ?feature_opp enpkg:has_retention_time ?rt_opp .
                ?feature_opp enpkg:has_sirius_annotation ?sirius_opp .
                    ?sirius_opp enpkg:has_InChIkey2D ?ik2d .

        FILTER(((?rt - 0.05) < ?rt_opp) && ((?rt + 0.05) > ?rt_opp))
    }
""")


results = sparql.queryAndConvert()
df = json_normalize(results['results']["bindings"])
df = df.stack().str.replace(uri_ref, "").unstack()
df.drop(list(df.filter(regex = 'type')), axis = 1, inplace = True)
df.columns = df.columns.str.replace('.value', '')
df
