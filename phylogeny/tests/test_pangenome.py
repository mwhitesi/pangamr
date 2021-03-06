
import os
import sys
import unittest

sys.path.insert(0, os.path.abspath('..'))

from panseqtrees import pangenome


class LocusTests(unittest.TestCase):

	def test_new_instance(self):

		l = pangenome.build_tree(
			{
				'first': 'ACAGCCATTGCCTGGTTGCTTCATGGGCAAAAGCTTGATGCGTGGGGCTTTGTAGGTATGGGGCTCATAATTGCTGCCTTTTTGCTCGCCCGATCCCCATCGTGGAAGTCGCTGCGGAGGCCGACGCCATGGTGACGGTGTTCGGCATTCTGAATCTCACCGAGGACTCCTTCTTCGATGAGAGCCGGCGGCTAGACCCCGCCGGCGCTGTCACCGCGGCGATCGAAATGCTGCGAGTCGGATCAGACGTCGTGGATGTCGGACCGGCCGCCAGCCATCCGGACGCGAGGCCTGTATCGCCGGCCGATGAGATCAGACGTATTGCGCCGCTCTTAGACGCCCTGTCCGATCAGATGCACCGTGTTTCAATCGACAGCTTCCAACCGGAAACCCAGCGCTATGCGCTCAAGCGCGGCGTGGGCTACCTGAACGATATCCAAGGATTTCCTGACCCTGCGCTCTATCCCGATATTGCTGAGGCGGACTGCAGGCTGGTGGTTATGCACTCAGCGCAGCGGGATGGCATCGCCACCCGCACCGGTCACCTTCGACCCGAAGACGCGCTCGACGAGATTGTGCGGTTCTTCGAGGCGCGGGTTTCCGCCTTGCGACGGAGCGGGGTCGCTGCCGACCGGCTCATCCTCGATCCGGGGATGGGATTTTTCTTGAGCCCCGCACCGGAAACATCGCTGCACGTGCTGTCGAACCTTCAAAAGCTGAAGTCGGCGTTGGGGCTTCCGCTATTGGTCTCGGTGTCGCGGAAATCCTTCTTGGGCGCCACCGTTGGCCTTCCTGTAAAGGATCTGGGTCCAGCGAGCCTTGCGGCGGAACTTCACGCGATCGGCAATGGCGCTGACTACGTCCGCACCCACGCGCCTGGAGATCTGCGAAGCGCAATCACCTTCTCGGAAACCCTCGCGAAATTTCGCAGTCGCGACGCCAGAGACCGAGGGTTAGATCATGCCTAGCATTC',
				'second': 'ACAGCCATTGCCTGGTTGCTTCATGGGCAAAAGCTTGATGCGTGGGGCTTTGTAGGTATGGGGCTCATAATTGCTGCCTTTTTGCTCGCCCGATCCCCATCGTGGAAGTCGCTGCGGAGGCCGACGCCATGGTGACGGTGTTCGGCATTCTGAATCTCACCGAGGACTCCTTCTTCGATGAGAGCCGGCGGCTAGACCCCGCCGGCGCTGTCACCGCGGCGATCGAAATGCTGCGAGTCGGATCAGACGTCGTGGATGTCGGACCGGCCGCCAGCCATCCGGACGCGAGGCCTGTATCGCCGGCCGATGAGATCAGACGTATTGCGCCGCTCTTAGACGCCCTGTCCGATCAGATGCACCGTGTTTCAATCGACAGCTTCCAACCGGAAACCCAGCGCTATGCGCTCAAGCGCGGCGTGGGCTACCTGAACGATATCCAAGGATTTCCTGACCCTGCGCTCTATCCCGATATTGCTGAGGCGGACTGCAGGCTGGTGGTTATGCACTCAGCGCAGCGGGATGGCATCGCCACCCGCACCGGTCACCTTCGACCCGAAGACGCGCTCGACGAGATTGTGCGGTTCTTCGAGGCGCGGGTTTCCGCCTTGCGACGGAGCGGGGTCGCTGCCGACCGGCTCATCCTCGATCCGGGGATGGGATTTTTCTTGAGCCCCGCACCGGAAACATCGCTGCACGTGCTGTCGAACCTTCAAAAGCTGAAGTCGGCGTTGGGGCTTCCGCTATTGGTCTCGGTGTCGCGGAAATCCTTCTTGGGCGCCACCGTTGGCCTTCCTGTAAAGGATCTGGGTCCAGCGAGCCTTGCGGCGGAACTTCACGCGATCGGCAATGGCGCTGACTACGTCCGCACCCACGCGCCTGGAGATCTGCGAAGCGCAATCACCTTCTCGGAAACCCTCGCGAAATTTCGCAGTCGCGACGCCAGAGACCGAGGGTTAGATCATGCCTAGCATTC',
				'third': 'ACAGCCATTGCCTGGTTGCTTCATGGGCAAAAGCTTGATGCGTGGGGCTTTGTAGGTATGGGGCTCATAATTGCTGCCTTTTTGCTCGCCCGATCCCCATCGTGGAAGTCGCTGCGGAGGCCGACGCCATGGTGACGGTGTTCGGCATTCTGAATCTCACCGAGGACTCCTTCTTCGATGAGAGCCGGCGGCTAGACCCCGCCGGCGCTGTCACCGCGGCGATCGAAATGCTGCGAGTCGGATCAGACGTCGTGGATGTCGGACCGGCCGCCAGCCATCCGGACGCGAGGCCTGTATCGCCGGCCGATGAGATCAGACGTATTGCGCCGCTCTTAGACGCCCTGTCCGATCAGATGCACCGTGTTTCAATCGACAGCTTCCAACCGGAAACCCAGCGCTATGCGCTCAAGCGCGGCGTGGGCTACCTGAACGATATCCAAGGATTTCCTGACCCTGCGCTCTATCCCGATATTGCTGAGGCGGACTGCAGGCTGGTGGTTATGCACTCAGCGCAGCGGGATGGCATCGCCACCCGCACCGGTCACCTTCGACCCGAAGACGCGCTCGACGAGATTGTGCGGTTCTTCGAGGCGCGGGTTTCCGCCTTGCGACGGAGCGGGGTCGCTGCCGACCGGCTCATCCTCGATCCGGGGATGGGATTTTTCTTGAGCCCCGCACCGGAAACATCGCTGCACGTGCTGTCGAACCTTCAAAAGCTGAAGTCGGCGTTGGGGCTTCCGCTATTGGTCTCGGTGTCGCGGAAATCCTTCTTGGGCGCCACCGTTGGCCTTCCTGTAAAGGATCTGGGTCCAGCGAGCCTTGCGGCGGAACTTCACGCGATCGGCAATGGCGCTGACTACGTCCGCACCCACGCGCCTGGAGATCTGCGAAGCGCAATCACCTTCTCGGAAACCCTCGCGAAATTTCGCAGTCGCGACGCCAGAGACCGAGGGTTAGATCATGCCTAGCATTC',
				'forth': 'TCACCACCGACTATTTGCAACAGTGCCACAGCCATTGCCTGGTTGCTTCATGGGCAAAAGCTTGATGCGTGGGGCTTTGTAGGTATGGGGCTCATAATTGCTGCCTTTTTGCTCGCCCGATCCCCATCGTGGAAGTCGCTGCGGAGGCCGACGCCATGGTGACGGTGTTCGGCATTCTGAATCTCACCGAGGACTCCTTCTTCGATGAGAGCCGGCGGCTAGACCCCGCCGGCGCTGTCACCGCGGCGATCGAAATGCTGCGAGTCGGATCAGACGTCGTGGATGTCGGACCGGCCGCCAGCCATCCGGACGCGAGGCCTGTATCGCCGGCCGATGAGATCAGACGTATTGCGCCGCTCTTAGACGCCCTGTCCGATCAGATGCACCGTGTTTCAATCGACAGCTTCCAACCGGAAACCCAGCGCTATGCGCTCAAGCGCGGCGTGGGCTACCTGAACGATATCCAAGGATTTCCTGACCCTGCGCTCTATCCCGATATTGCTGAGGCGGACTGCAGGCTGGTGGTTATGCACTCAGCGCAGCGGGATGGCATCGCCACCCGCACCGGTCACCTTCGACCCGAAGACGCGCTCGACGAGATTGTGCGGTTCTTCGAGGCGCGGGTTTCCGCCTTGCGACGGAGCGGGGTCGCTGCCGACCGGCTCATCCTCGATCCGGGGATGGGATTTTTCTTGAGCCCCGCACCGGAAACATCGCTGCACGTGCTGTCGAACCTTCAAAAGCTGAAGTCGGCGTTGGGGCTTCCGCTATTGGTCTCGGTGTCGCGGAAATCCTTCTTGGGCGCCACCGTTGGCCTTCCTGTAAAGGATCTGGGTCCAGCGAGCCTTGCGGCGGAACTTCACGCGATCGGCAATGGCGCTGACTACGTCCGCACCCACGCGCCTGGAGATCTGCGAAGCGCAATCACCTTCTCGGAAACCCTCGCGAAATTTCGCAGTCGCGACGCCAGAGACCGAGGGTTAGATCATGCCTAGCATTC'
			}	
		)

		self.assertIsNotNone(l)
		self.assertTrue('second' in l._id2node.keys())
		self.assertTrue(l._id2node['second'].id == 'second')
		self.assertIsNotNone(l._tree)


	def test_new_instance_with_phenotype(self):

		l = pangenome.build_tree(
			{
				'first': 'ACAGCCATTGCCTGGTTGCTTCATGGGCAAAAGCTTGATGCGTGGGGCTTTGTAGGTATGGGGCTCATAATTGCTGCCTTTTTGCTCGCCCGATCCCCATCGTGGAAGTCGCTGCGGAGGCCGACGCCATGGTGACGGTGTTCGGCATTCTGAATCTCACCGAGGACTCCTTCTTCGATGAGAGCCGGCGGCTAGACCCCGCCGGCGCTGTCACCGCGGCGATCGAAATGCTGCGAGTCGGATCAGACGTCGTGGATGTCGGACCGGCCGCCAGCCATCCGGACGCGAGGCCTGTATCGCCGGCCGATGAGATCAGACGTATTGCGCCGCTCTTAGACGCCCTGTCCGATCAGATGCACCGTGTTTCAATCGACAGCTTCCAACCGGAAACCCAGCGCTATGCGCTCAAGCGCGGCGTGGGCTACCTGAACGATATCCAAGGATTTCCTGACCCTGCGCTCTATCCCGATATTGCTGAGGCGGACTGCAGGCTGGTGGTTATGCACTCAGCGCAGCGGGATGGCATCGCCACCCGCACCGGTCACCTTCGACCCGAAGACGCGCTCGACGAGATTGTGCGGTTCTTCGAGGCGCGGGTTTCCGCCTTGCGACGGAGCGGGGTCGCTGCCGACCGGCTCATCCTCGATCCGGGGATGGGATTTTTCTTGAGCCCCGCACCGGAAACATCGCTGCACGTGCTGTCGAACCTTCAAAAGCTGAAGTCGGCGTTGGGGCTTCCGCTATTGGTCTCGGTGTCGCGGAAATCCTTCTTGGGCGCCACCGTTGGCCTTCCTGTAAAGGATCTGGGTCCAGCGAGCCTTGCGGCGGAACTTCACGCGATCGGCAATGGCGCTGACTACGTCCGCACCCACGCGCCTGGAGATCTGCGAAGCGCAATCACCTTCTCGGAAACCCTCGCGAAATTTCGCAGTCGCGACGCCAGAGACCGAGGGTTAGATCATGCCTAGCATTC',
				'second': 'ACAGCCATTGCCTGGTTGCTTCATGGGCAAAAGCTTGATGCGTGGGGCTTTGTAGGTATGGGGCTCATAATTGCTGCCTTTTTGCTCGCCCGATCCCCATCGTGGAAGTCGCTGCGGAGGCCGACGCCATGGTGACGGTGTTCGGCATTCTGAATCTCACCGAGGACTCCTTCTTCGATGAGAGCCGGCGGCTAGACCCCGCCGGCGCTGTCACCGCGGCGATCGAAATGCTGCGAGTCGGATCAGACGTCGTGGATGTCGGACCGGCCGCCAGCCATCCGGACGCGAGGCCTGTATCGCCGGCCGATGAGATCAGACGTATTGCGCCGCTCTTAGACGCCCTGTCCGATCAGATGCACCGTGTTTCAATCGACAGCTTCCAACCGGAAACCCAGCGCTATGCGCTCAAGCGCGGCGTGGGCTACCTGAACGATATCCAAGGATTTCCTGACCCTGCGCTCTATCCCGATATTGCTGAGGCGGACTGCAGGCTGGTGGTTATGCACTCAGCGCAGCGGGATGGCATCGCCACCCGCACCGGTCACCTTCGACCCGAAGACGCGCTCGACGAGATTGTGCGGTTCTTCGAGGCGCGGGTTTCCGCCTTGCGACGGAGCGGGGTCGCTGCCGACCGGCTCATCCTCGATCCGGGGATGGGATTTTTCTTGAGCCCCGCACCGGAAACATCGCTGCACGTGCTGTCGAACCTTCAAAAGCTGAAGTCGGCGTTGGGGCTTCCGCTATTGGTCTCGGTGTCGCGGAAATCCTTCTTGGGCGCCACCGTTGGCCTTCCTGTAAAGGATCTGGGTCCAGCGAGCCTTGCGGCGGAACTTCACGCGATCGGCAATGGCGCTGACTACGTCCGCACCCACGCGCCTGGAGATCTGCGAAGCGCAATCACCTTCTCGGAAACCCTCGCGAAATTTCGCAGTCGCGACGCCAGAGACCGAGGGTTAGATCATGCCTAGCATTC',
				'third': 'ACAGCCATTGCCTGGTTGCTTCATGGGCAAAAGCTTGATGCGTGGGGCTTTGTAGGTATGGGGCTCATAATTGCTGCCTTTTTGCTCGCCCGATCCCCATCGTGGAAGTCGCTGCGGAGGCCGACGCCATGGTGACGGTGTTCGGCATTCTGAATCTCACCGAGGACTCCTTCTTCGATGAGAGCCGGCGGCTAGACCCCGCCGGCGCTGTCACCGCGGCGATCGAAATGCTGCGAGTCGGATCAGACGTCGTGGATGTCGGACCGGCCGCCAGCCATCCGGACGCGAGGCCTGTATCGCCGGCCGATGAGATCAGACGTATTGCGCCGCTCTTAGACGCCCTGTCCGATCAGATGCACCGTGTTTCAATCGACAGCTTCCAACCGGAAACCCAGCGCTATGCGCTCAAGCGCGGCGTGGGCTACCTGAACGATATCCAAGGATTTCCTGACCCTGCGCTCTATCCCGATATTGCTGAGGCGGACTGCAGGCTGGTGGTTATGCACTCAGCGCAGCGGGATGGCATCGCCACCCGCACCGGTCACCTTCGACCCGAAGACGCGCTCGACGAGATTGTGCGGTTCTTCGAGGCGCGGGTTTCCGCCTTGCGACGGAGCGGGGTCGCTGCCGACCGGCTCATCCTCGATCCGGGGATGGGATTTTTCTTGAGCCCCGCACCGGAAACATCGCTGCACGTGCTGTCGAACCTTCAAAAGCTGAAGTCGGCGTTGGGGCTTCCGCTATTGGTCTCGGTGTCGCGGAAATCCTTCTTGGGCGCCACCGTTGGCCTTCCTGTAAAGGATCTGGGTCCAGCGAGCCTTGCGGCGGAACTTCACGCGATCGGCAATGGCGCTGACTACGTCCGCACCCACGCGCCTGGAGATCTGCGAAGCGCAATCACCTTCTCGGAAACCCTCGCGAAATTTCGCAGTCGCGACGCCAGAGACCGAGGGTTAGATCATGCCTAGCATTC',
				'forth': 'TCACCACCGACTATTTGCAACAGTGCCACAGCCATTGCCTGGTTGCTTCATGGGCAAAAGCTTGATGCGTGGGGCTTTGTAGGTATGGGGCTCATAATTGCTGCCTTTTTGCTCGCCCGATCCCCATCGTGGAAGTCGCTGCGGAGGCCGACGCCATGGTGACGGTGTTCGGCATTCTGAATCTCACCGAGGACTCCTTCTTCGATGAGAGCCGGCGGCTAGACCCCGCCGGCGCTGTCACCGCGGCGATCGAAATGCTGCGAGTCGGATCAGACGTCGTGGATGTCGGACCGGCCGCCAGCCATCCGGACGCGAGGCCTGTATCGCCGGCCGATGAGATCAGACGTATTGCGCCGCTCTTAGACGCCCTGTCCGATCAGATGCACCGTGTTTCAATCGACAGCTTCCAACCGGAAACCCAGCGCTATGCGCTCAAGCGCGGCGTGGGCTACCTGAACGATATCCAAGGATTTCCTGACCCTGCGCTCTATCCCGATATTGCTGAGGCGGACTGCAGGCTGGTGGTTATGCACTCAGCGCAGCGGGATGGCATCGCCACCCGCACCGGTCACCTTCGACCCGAAGACGCGCTCGACGAGATTGTGCGGTTCTTCGAGGCGCGGGTTTCCGCCTTGCGACGGAGCGGGGTCGCTGCCGACCGGCTCATCCTCGATCCGGGGATGGGATTTTTCTTGAGCCCCGCACCGGAAACATCGCTGCACGTGCTGTCGAACCTTCAAAAGCTGAAGTCGGCGTTGGGGCTTCCGCTATTGGTCTCGGTGTCGCGGAAATCCTTCTTGGGCGCCACCGTTGGCCTTCCTGTAAAGGATCTGGGTCCAGCGAGCCTTGCGGCGGAACTTCACGCGATCGGCAATGGCGCTGACTACGTCCGCACCCACGCGCCTGGAGATCTGCGAAGCGCAATCACCTTCTCGGAAACCCTCGCGAAATTTCGCAGTCGCGACGCCAGAGACCGAGGGTTAGATCATGCCTAGCATTC'
			},
			phenodict = {
				'pheno1': {
					'first': 'Resistant',
					'second': 'Susceptible',
					'third': 'Resistant'
				}
			}
		)

		self.assertIsNotNone(l)
		self.assertIsNotNone(l._pheno_classes)
		self.assertTrue(len(l._pheno_classes)==1)
		self.assertTrue(len(l._pheno_classes['pheno1'])==2)
		self.assertEqual(l._id2node['second'].pheno1, 'Susceptible')
		leafnode = l._tree.get_leaves_by_name(l._id2node['forth'].name)[0]
		self.assertEqual(leafnode.pheno1, l._pheno_na)


	def test_ncolors(self):

		colors = pangenome.ncolors([2,1])

		self.assertTrue(len(colors) == 2)
		self.assertTrue(len(colors[0]) == 2)
		self.assertTrue(len(colors[1]) == 1)


def main():
    unittest.main()

if __name__ == '__main__':
    main()
	