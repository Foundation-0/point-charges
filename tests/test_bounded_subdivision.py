import importlib.util,pathlib,sys,unittest,numpy as np
ROOT=pathlib.Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT));import certify
class BoundedSubdivision(unittest.TestCase):
 def test_small_batches_preserve_three_charge_count(self):
  sites=np.diag([3.,3.,3.]);charges=np.ones(3)
  for batch in (8,32,8192):
   zeros,stats=certify.certify(sites,charges,batch_size=batch,verbose=False)
   self.assertEqual(len(zeros),4);self.assertEqual(sorted(z['index']for z in zeros),[1,2,2,2]);self.assertGreater(stats['boxes'],0)
 def test_invalid_batch(self):
  for batch in (0,-1,True):
   with self.assertRaises(ValueError):certify.certify(np.diag([3.,3.,3.]),np.ones(3),batch_size=batch,verbose=False)
 def test_budget_does_not_emit_result(self):
  with self.assertRaisesRegex(RuntimeError,'box budget'):certify.certify(np.diag([3.,3.,3.]),np.ones(3),batch_size=8,max_boxes=1,verbose=False)
if __name__=='__main__':unittest.main()
