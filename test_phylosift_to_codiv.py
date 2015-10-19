#!/usr/bin/env python

from StringIO import StringIO
from unittest import TestCase, main
import pandas as pd
import pandas.util.testing as pdt

from phylosift_to_codiv import parse_ps_taxa_summary, parse_lencovgc


class PhylosiftTests(TestCase):
    def setUp(self):
        print "starting test"

    def test_parse_ps_taxa_summary(self):
        obs_taxa_df = parse_ps_taxa_summary(StringIO(taxa_summary_file))

        pdt.assert_frame_equal(obs_taxa_df, pd.DataFrame.from_dict(exp_taxa_dict, orient='index'))

    def test_parse_lencovgc(self):
        obs_lencovgc_df = parse_lencovgc(StringIO(lencovgc_file))

        pdt.assert_frame_equal(obs_lencovgc_df, exp_lencovgc_df)



exp_taxa_dict = {
'scaffold_124805_325_3.3754_0.5477': {'superkingdom': 'BACTERIA', 'phylum': 'PROTEOBACTERIA', 'class': 'ALPHAPROTEOBACTERIA', 'order': 'RHIZOBIALES', 'family': 'BRUCELLACEAE', 'genus': 'BRUCELLA', 'species': None},
'scaffold_10482_1294_7.3223_0.6561': {'superkingdom': 'BACTERIA', 'phylum': 'PROTEOBACTERIA', 'class': 'BETAPROTEOBACTERIA', 'order': 'BURKHOLDERIALES', 'family': 'ALCALIGENACEAE', 'genus': None, 'species': None},
'scaffold_63179_551_2.1561_0.4120': {'superkingdom': 'EUKARYOTA', 'phylum': None, 'class': None, 'order': None, 'family': None, 'genus': None, 'species': None},
'scaffold_3006_2632_15.6269_0.6045': {'superkingdom': 'BACTERIA', 'phylum': 'PROTEOBACTERIA', 'class': 'BETAPROTEOBACTERIA', 'order': 'BURKHOLDERIALES', 'family': 'ALCALIGENACEAE', 'genus': None, 'species': None},
'scaffold_81902_485_4.0701_0.6124': {'superkingdom': 'ARCHAEA', 'phylum': 'THAUMARCHAEOTA', 'class': None, 'order': 'NITROSOPUMILALES', 'family': 'NITROSOPUMILACEAE', 'genus': 'NITROSOPUMILUS', 'species': None}
}

lencovgc_file = """JS0206_persimilis	scaffold_124805	325	3.37538461538462	0.547692307692308
JS0206_persimilis	scaffold_10482	1294	7.32225656877899	0.656105100463679
JS0206_persimilis	scaffold_63179	551	2.15607985480944	0.411978221415608
JS0206_persimilis	scaffold_3006	2632	15.6268996960486	0.604483282674772
JS0206_persimilis	scaffold_81902	485	4.07010309278351	0.612371134020619"""

exp_lencovgc_df = pd.DataFrame({'len': [1294, 325, 2632, 551, 485],
'gc': [0.656105100463679, 0.547692307692308, 0.604483282674772, 0.411978221415608, 0.612371134020619],
'seq': ['scaffold_10482', 'scaffold_124805', 'scaffold_3006', 'scaffold_63179', 'scaffold_81902'],
'file': ['JS0206_persimilis', 'JS0206_persimilis', 'JS0206_persimilis', 'JS0206_persimilis', 'JS0206_persimilis'],
'cov': [7.32225656877899, 3.37538461538462, 15.6268996960486, 2.15607985480944, 4.07010309278351]},
index = ['scaffold_10482_1294_7.3223_0.6561',
         'scaffold_124805_325_3.3754_0.5477',
         'scaffold_3006_2632_15.6269_0.6045',
         'scaffold_63179_551_2.1561_0.4120',
         'scaffold_81902_485_4.0701_0.6124'],
columns = ['file','seq','len','cov','gc']).reindex(
        ['scaffold_124805_325_3.3754_0.5477',
         'scaffold_10482_1294_7.3223_0.6561',
         'scaffold_63179_551_2.1561_0.4120',
         'scaffold_3006_2632_15.6269_0.6045',
         'scaffold_81902_485_4.0701_0.6124'])

taxa_summary_file = """#Sequence_ID	Hit_Coordinates	NCBI_Taxon_ID	Taxon_Rank	Taxon_Name	Cumulative_Probability_Mass	Markers_Hit
scaffold_124805_325_3.3754_0.5477	1.270	1224	phylum	PROTEOBACTERIA	1	concat
scaffold_124805_325_3.3754_0.5477	1.270	234	genus	BRUCELLA	1	concat
scaffold_124805_325_3.3754_0.5477	1.270	28211	class	ALPHAPROTEOBACTERIA	1	concat
scaffold_124805_325_3.3754_0.5477	1.270	2	superkingdom	BACTERIA	1	concat
scaffold_124805_325_3.3754_0.5477	1.270	356	order	RHIZOBIALES	1	concat
scaffold_124805_325_3.3754_0.5477	1.270	1	no rank	ROOT	1	concat
scaffold_124805_325_3.3754_0.5477	1.270	118882	family	BRUCELLACEAE	1	concat
scaffold_124805_325_3.3754_0.5477	1.270	131567	no rank	CELLULAR ORGANISMS	1	concat
scaffold_124805_325_3.3754_0.5477	1.270	236	species	BRUCELLA OVIS	0.146324604707971	concat
scaffold_124805_325_3.3754_0.5477	1.270	1218315	species	BRUCELLA INOPINATA	0.122547863410975	concat
scaffold_124805_325_3.3754_0.5477	1.270	693750	species	BRUCELLA SP. BO2	0.0993315903380801	concat
scaffold_124805_325_3.3754_0.5477	1.270	444178	no rank	BRUCELLA OVIS ATCC 25840	0.0731623023539854	concat
scaffold_124805_325_3.3754_0.5477	1.270	470735	no rank	BRUCELLA INOPINATA BO1	0.0612739317054877	concat
scaffold_124805_325_3.3754_0.5477	1.270	693748	species	BRUCELLA SP. NF 2653	0.0409313561130918	concat
scaffold_10482_1294_7.3223_0.6561	3.1214	1	no rank	ROOT	1	concat
scaffold_10482_1294_7.3223_0.6561	3.1214	506	family	ALCALIGENACEAE	1	concat
scaffold_10482_1294_7.3223_0.6561	3.1214	1224	phylum	PROTEOBACTERIA	1	concat
scaffold_10482_1294_7.3223_0.6561	3.1214	131567	no rank	CELLULAR ORGANISMS	1	concat
scaffold_10482_1294_7.3223_0.6561	3.1214	2	superkingdom	BACTERIA	1	concat
scaffold_10482_1294_7.3223_0.6561	3.1214	28216	class	BETAPROTEOBACTERIA	1	concat
scaffold_10482_1294_7.3223_0.6561	3.1214	80840	order	BURKHOLDERIALES	1	concat
scaffold_63179_551_2.1561_0.4120	318.551	1	no rank	ROOT	1	concat
scaffold_63179_551_2.1561_0.4120	318.551	131567	no rank	CELLULAR ORGANISMS	1	concat
scaffold_63179_551_2.1561_0.4120	318.551	2759	superkingdom	EUKARYOTA	0.914458026858309	concat
scaffold_63179_551_2.1561_0.4120	318.551	33154	no rank	OPISTHOKONTA	0.809280286576919	concat
scaffold_63179_551_2.1561_0.4120	318.551	4751	kingdom	FUNGI	0.809280286576919	concat
scaffold_63179_551_2.1561_0.4120	318.551	33634	no rank	STRAMENOPILES	0.105177740281389	concat
scaffold_63179_551_2.1561_0.4120	318.551	2	superkingdom	BACTERIA	0.0855419731416915	concat
scaffold_3006_2632_15.6269_0.6045	396.758.2359.2631.1.345	1	no rank	ROOT	1	concat
scaffold_3006_2632_15.6269_0.6045	396.758.2359.2631.1.345	506	family	ALCALIGENACEAE	1	concat
scaffold_3006_2632_15.6269_0.6045	396.758.2359.2631.1.345	1224	phylum	PROTEOBACTERIA	1	concat
scaffold_3006_2632_15.6269_0.6045	396.758.2359.2631.1.345	131567	no rank	CELLULAR ORGANISMS	1	concat
scaffold_3006_2632_15.6269_0.6045	396.758.2359.2631.1.345	2	superkingdom	BACTERIA	1	concat
scaffold_3006_2632_15.6269_0.6045	396.758.2359.2631.1.345	28216	class	BETAPROTEOBACTERIA	1	concat
scaffold_3006_2632_15.6269_0.6045	396.758.2359.2631.1.345	80840	order	BURKHOLDERIALES	1	concat
scaffold_81902_485_4.0701_0.6124	96.476	2157	superkingdom	ARCHAEA	1	concat
scaffold_81902_485_4.0701_0.6124	96.476	1	no rank	ROOT	1	concat
scaffold_81902_485_4.0701_0.6124	96.476	31932	order	NITROSOPUMILALES	1	concat
scaffold_81902_485_4.0701_0.6124	96.476	131567	no rank	CELLULAR ORGANISMS	1	concat
scaffold_81902_485_4.0701_0.6124	96.476	651137	phylum	THAUMARCHAEOTA	1	concat
scaffold_81902_485_4.0701_0.6124	96.476	338190	family	NITROSOPUMILACEAE	1	concat
scaffold_81902_485_4.0701_0.6124	96.476	338191	genus	NITROSOPUMILUS	0.67901729714	concat
scaffold_81902_485_4.0701_0.6124	96.476	1007084	species	CANDIDATUS NITROSOARCHAEUM LIMNIA	0.32098270286	concat
scaffold_81902_485_4.0701_0.6124	96.476	1007082	genus	CANDIDATUS NITROSOARCHAEUM	0.32098270286	concat
scaffold_81902_485_4.0701_0.6124	96.476	1027373	species	NITROSOPUMILUS SP. AR	0.171759969985	concat
scaffold_81902_485_4.0701_0.6124	96.476	1027374	species	NITROSOPUMILUS SP. SJ	0.171759969985	concat
scaffold_81902_485_4.0701_0.6124	96.476	859192	no rank	CANDIDATUS NITROSOARCHAEUM LIMNIA BG20	0.16049135143	concat
scaffold_81902_485_4.0701_0.6124	96.476	886738	no rank	CANDIDATUS NITROSOARCHAEUM LIMNIA SFB1	0.16049135143	concat
scaffold_81902_485_4.0701_0.6124	96.476	338192	species	NITROSOPUMILUS MARITIMUS	0.116313762662	concat
scaffold_81902_485_4.0701_0.6124	96.476	1229909	species	CANDIDATUS NITROSOPUMILUS SP. AR2	0.102869831846	concat
scaffold_81902_485_4.0701_0.6124	96.476	436308	no rank	NITROSOPUMILUS MARITIMUS SCM1	0.058156881331	concat
"""


# run tests if called from command line
if __name__ == "__main__":
    main()
