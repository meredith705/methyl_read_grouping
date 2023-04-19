ONT Methylation Read Clustering

Motivation: 
CARD brain samples are compromised of several cell types. Methylation atlas of pure cell types paper prove that there are regions of specific hypo/hyper methylation per cell type ( https://www.nature.com/articles/s41586-022-05580-6 ). 
Proof of Principle paper: 
https://www.nature.com/articles/s41588-022-01188-8


Aims
This work aims to identify cell types, or at least the number of potential cell types, within a sample based off the methylation profile of each read. By identifying reads that come from particular cell types can allow us to have a more fine grained analysis of differential methylation within each cell type in response to control/pathogenic brain samples. 
Note: diseased brains may complicate this analysis as their methylation profiles may be altered by disease. 

Methods: 
Parse the bam and get all the methyl sites in terms of reference cooridnates. 
Build an overlap graph of the CpG sites 
	Calculate overlap with a model that accounts for
		Measurement error 
		Information at each site ( degree of mixed methylation at each stie) 
		Remove isolated vertices 
		Keep significant overlaps constructs local ‘methylotpyes’ 
Use this graph to infer cell types within the sample 





63678769 11 dict_keys([
'4cb0c950-c48c-4fc9-b383-7af464e2c040', +
'fbd13903-485c-4fed-89a8-4dc1975ebd71', +
'4c5c8a96-c1f7-4220-a423-0071d80e809c', +
'e8bec659-5b84-45b8-8116-a909a08d05e2', +
'52bd741f-90f4-43f0-af37-d4d0b6159ed7', +
'a523e12c-b6f4-46a0-ae33-31d58eecee8a', +
'500c6b22-99e1-4400-b532-883480152cba', +
'2283bc65-9368-4105-9e85-74fee5905a85', +
'cfc5bdb5-577a-459f-a329-1d262dc6ef53', +
'30e1d055-82ff-4ebc-9773-47076f6cab0d', +
'a90cd397-9b2b-499f-8c19-98099c4a362d']) +

63678771 7 dict_keys([
'7bb9b12c-08b5-44f2-bd2b-8e37fbc27eff',  -
'863bf4ab-2057-4f93-a86b-946dd3f08331', -
'32fb2648-7f50-4b11-8f56-7a43539901ab', -
'2e42150e-1d8d-4039-bbf7-24df2fc74a93', -
'fc5046ce-f08d-498a-9bfe-583bdbb8f338', -
'4e522e43-10c2-41d5-97a7-a8711ac7c0e7', -
'7e874550-947e-48b7-8562-e922645a70c5']) -


wider fetch batch 
63678769 28 dict_keys(['4cb0c950-c48c-4fc9-b383-7af464e2c040', 'fbd13903-485c-4fed-89a8-4dc1975ebd71', '4c5c8a96-c1f7-4220-a423-0071d80e809c', 'e8bec659-5b84-45b8-8116-a909a08d05e2', '52bd741f-90f4-43f0-af37-d4d0b6159ed7', 'a523e12c-b6f4-46a0-ae33-31d58eecee8a', '500c6b22-99e1-4400-b532-883480152cba', '2283bc65-9368-4105-9e85-74fee5905a85', 'cfc5bdb5-577a-459f-a329-1d262dc6ef53', '30e1d055-82ff-4ebc-9773-47076f6cab0d', 'a90cd397-9b2b-499f-8c19-98099c4a362d', '9d57463c-fe29-4b1e-805c-0b2842d2fb07', '87050bca-dc81-43f8-8149-943e151aa6f9', 'fbb3fcd8-d470-468e-801d-ccd4889a67e5', 'ab5d7f5b-bbd4-45ef-8d26-e992430259f3', 'bcb241f5-1e58-42ff-af0c-cf051ebb6f1a', '68cdff21-3dbc-4434-9d23-0c7bc5ba3ef7', 'a1b6d931-0a5a-4a0d-8ba5-70cc238ff0c6', '0060f390-b24e-4796-ab33-014022fbd8fd', '6fff6026-fcbb-4716-a55c-4309cf237d20', '270d2664-95e1-464c-b6b2-df30f18a53d8', '5e398f24-902c-4151-8830-815e583933eb', 'ff2ff66a-32b5-4e78-a2f8-b5f655aabbf6', 'b475fc05-9534-4530-a972-4515ae8c1449', 'e60be7d4-c890-404c-bff9-fb90b5e767c4', '75e87014-9e99-46c1-aec0-67da4ca4bd78', '9f0afe59-004c-44d6-bda4-039270b25ec6', 
'c26681c3-18d0-4695-97c1-a6a24426b8d1'])

63678771 18 dict_keys(['7bb9b12c-08b5-44f2-bd2b-8e37fbc27eff', '863bf4ab-2057-4f93-a86b-946dd3f08331', '32fb2648-7f50-4b11-8f56-7a43539901ab', '2e42150e-1d8d-4039-bbf7-24df2fc74a93', 'fc5046ce-f08d-498a-9bfe-583bdbb8f338', '4e522e43-10c2-41d5-97a7-a8711ac7c0e7', '7e874550-947e-48b7-8562-e922645a70c5', '99a8fff2-a630-426d-92ed-1b567cb92cc6', '1345a9f8-e8b4-4c12-95f1-0496c1aa625a', 'caa7ee63-3e63-4547-810b-6d33c3b70b72', '891c5ed6-705b-4bb1-b0c3-fc9c9aa24817', '86024447-1b5d-4827-ba90-564fbc50dee0', '6cfa12d4-b109-4aca-a2aa-32b29efc7ea9', '0bedcf21-61b1-4cc7-80f6-1436b2c14222', '2fb2b4dd-6984-48d8-81dc-8bbaee8bfe12', 'f3f4161b-88af-4840-9e8d-30b563ee6952', 'bb1b8bff-e203-457b-babf-ba7211305325', '192759ef-b07b-4be2-869a-d62c034891a7'])
