import unittest
import numpy as np
import mpT3libs

class mpT3libsUnitTests(unittest.TestCase):

	# OUTPUT FINITENESS TESTS
	def test_FiniteValues_H3_rnd(self):
		"""
		Testing if H3_rnd(shp3) returns FINITE values (no nan/inf allowed)
		"""
		maxshp = 16
		shp3 = np.array(np.random.random_integers(maxshp,high=None,size=3))
		
		H11,H12,H13,H22,H23,H33 = mpT3libs.H3_rnd( shp3 )
		
		self.assertTrue( np.all( np.isfinite( H11 ) ) &
						 np.all( np.isfinite( H12 ) ) &
						 np.all( np.isfinite( H13 ) ) &
						 np.all( np.isfinite( H22 ) ) &
						 np.all( np.isfinite( H23 ) ) &
						 np.all( np.isfinite( H33 ) ) )

	def test_FiniteValues_T3_iso(self):
		"""
		Testing if T3_iso(shp3) returns FINITE values (no nan/inf allowed)
		"""
		maxshp = 16
		shp3 = np.array(np.random.random_integers(maxshp,high=None,size=3))

		El1,El2,El3,Ev1,Ev2,Ev3 = mpT3libs.T3_iso( shp3 )

		self.assertTrue( np.all( np.isfinite( El1 ) ) &
						 np.all( np.isfinite( El2 ) ) &
						 np.all( np.isfinite( El3 ) ) &
						 np.all( np.isfinite( Ev1 ) ) &
						 np.all( np.isfinite( Ev2 ) ) &
						 np.all( np.isfinite( Ev3 ) ) )

	def test_FiniteValues_T3_rnd(self):
		"""
		Testing if T3_rnd(shp3) returns FINITE values (no nan/inf allowed)
		"""
		maxshp = 16
		shp3 = np.array(np.random.random_integers(maxshp,high=None,size=3))

		El1,El2,El3,Ev1,Ev2,Ev3 = mpT3libs.T3_rnd( shp3 )

		self.assertTrue( np.all( np.isfinite( El1 ) ) &
						 np.all( np.isfinite( El2 ) ) &
						 np.all( np.isfinite( El3 ) ) &
						 np.all( np.isfinite( Ev1 ) ) &
						 np.all( np.isfinite( Ev2 ) ) &
						 np.all( np.isfinite( Ev3 ) ) )

	def test_FiniteValues_msk3_to_idx_full(self):
		"""
		Testing if msk3_to_idx(msk3) returns FINITE values (no nan/inf allowed)
		with msk3 a boolean array identically True
		"""
		maxshp = 16
		shp3 = np.array(np.random.random_integers(maxshp,high=None,size=3))

		msk3_full = np.ones( shp3 , dtype=bool )

		idx_full = mpT3libs.msk3_to_idx( msk3_full )

		self.assertTrue( np.all( np.isfinite( idx_full ) ) &
						 np.all( np.prod(np.shape(idx_full)) == np.sum(msk3_full) ) )

	def test_FiniteValues_msk3_to_idx_part(self):
		"""
		Testing if msk3_to_idx(msk3) returns FINITE values (no nan/inf allowed)
		with msk3 a boolean array with random True and False values
		"""
		maxshp = 16
		shp3 = np.array(np.random.random_integers(maxshp,high=None,size=3))

		msk3_part = np.random.rand( shp3[0], shp3[1], shp3[2] ) < 0.5

		idx_part = mpT3libs.msk3_to_idx( msk3_part )

		self.assertTrue( np.all( np.isfinite( idx_part ) ) &
						 np.all( np.prod(np.shape(idx_part)) == np.sum(msk3_part) ) )

	def test_FiniteValues_msk3_to_idx_null(self):
		"""
		Testing if msk3_to_idx(msk3) returns FINITE values (no nan/inf allowed)
		vith msk3 a boolean array identically False
		"""
		maxshp = 16
		shp3 = np.array(np.random.random_integers(maxshp,high=None,size=3))

		msk3_null = np.zeros( shp3 , dtype=bool )

		idx_null = mpT3libs.msk3_to_idx( msk3_null )

		self.assertTrue( np.all( np.isfinite( idx_null ) ) &
						 np.all( np.prod(np.shape(idx_null)) == np.sum(msk3_null) ) )

	def test_FiniteValues_iniElvs(self):
		"""
		Testing if iniElvs(shp3) returns FINITE values (no nan/inf allowed)
		"""
		maxshp = 16
		shp3 = np.array(np.random.random_integers(maxshp,high=None,size=3))

		(El1,El2,El3,
		 Ev11,Ev12,Ev13,
		 Ev21,Ev22,Ev23,
		 Ev31,Ev32,Ev33) = mpT3libs.iniElvs( shp3 )

		self.assertTrue( np.all( np.isfinite( El1 ) ) &
						 np.all( np.isfinite( El2 ) ) &
						 np.all( np.isfinite( El3 ) ) &
						 np.all( np.isfinite( Ev11 ) ) &
						 np.all( np.isfinite( Ev12 ) ) &
						 np.all( np.isfinite( Ev13 ) ) &
						 np.all( np.isfinite( Ev21 ) ) &
						 np.all( np.isfinite( Ev22 ) ) &
						 np.all( np.isfinite( Ev23 ) ) &
						 np.all( np.isfinite( Ev31 ) ) &
						 np.all( np.isfinite( Ev32 ) ) &
						 np.all( np.isfinite( Ev33 ) ) )

	def test_FiniteValues_catEvs(self):
		"""
		Testing if catEvs(Ev11,Ev12,...,Ev32,Ev33) returns FINITE values (no nan/inf allowed)
		Testing if returned Ev1,Ev2,Ev3 are 4D
		"""
		maxshp = 16
		shp3 = np.array(np.random.random_integers(maxshp,high=None,size=3))

		( _ , _ , _ ,
		 Ev11,Ev12,Ev13,
		 Ev21,Ev22,Ev23,
		 Ev31,Ev32,Ev33) = mpT3libs.iniElvs( shp3 )

		Ev1,Ev2,Ev3 = mpT3libs.catEvs(Ev11,Ev12,Ev13,Ev21,Ev22,Ev23,Ev31,Ev32,Ev33)

		self.assertTrue( np.all( np.isfinite( Ev1 ) ) &
						 np.all( np.isfinite( Ev2 ) ) &
						 np.all( np.isfinite( Ev3 ) ) &
						 np.all( np.prod(np.shape(np.shape(Ev1))) == 4 ) &
						 np.all( np.prod(np.shape(np.shape(Ev2))) == 4 ) &
						 np.all( np.prod(np.shape(np.shape(Ev3))) == 4 ) )

	def test_FiniteValues_sepEvs(self):
		"""
		Testing if sepEvs(Ev1,Ev2,Ev3) returns FINITE values (no nan/inf allowed)
		Testing if returned Ev11,Ev12,...,Ev32,Ev33 are 3D
		"""
		maxshp = 16
		shp3 = np.array(np.random.random_integers(maxshp,high=None,size=3))

		_ , _ , _ ,Ev1,Ev2,Ev3 = mpT3libs.T3_rnd( shp3 )

		Ev11,Ev12,Ev13,Ev21,Ev22,Ev23,Ev31,Ev32,Ev33 = mpT3libs.sepEvs( Ev1,Ev2,Ev3 )

		self.assertTrue( np.all( np.isfinite( Ev11 ) ) &
						 np.all( np.isfinite( Ev12 ) ) &
						 np.all( np.isfinite( Ev13 ) ) &
						 np.all( np.isfinite( Ev21 ) ) &
						 np.all( np.isfinite( Ev22 ) ) &
						 np.all( np.isfinite( Ev23 ) ) &
						 np.all( np.isfinite( Ev31 ) ) &
						 np.all( np.isfinite( Ev32 ) ) &
						 np.all( np.isfinite( Ev33 ) ) &
						 np.all(np.shape(Ev11) == shp3) &
						 np.all(np.shape(Ev12) == shp3) &
						 np.all(np.shape(Ev13) == shp3) &
						 np.all(np.shape(Ev21) == shp3) &
						 np.all(np.shape(Ev22) == shp3) &
						 np.all(np.shape(Ev23) == shp3) &
						 np.all(np.shape(Ev31) == shp3) &
						 np.all(np.shape(Ev32) == shp3) &
						 np.all(np.shape(Ev33) == shp3) )

	def test_FiniteValues_rglElvs(self):
		"""
		Testing if rglElvs(El1,El2,El3,Ev1,Ev2,Ev3) returns FINITE values (no nan/inf allowed)
		Testing if type is np.single (float32)
		"""
		maxshp = 16
		shp3 = np.array(np.random.random_integers(maxshp,high=None,size=3))

		El1,El2,El3,Ev1,Ev2,Ev3 = mpT3libs.T3_rnd( shp3 )

		El1,El2,El3,Ev1,Ev2,Ev3 = mpT3libs.rglElvs(El1,El2,El3,Ev1,Ev2,Ev3)

		self.assertTrue( np.all( np.isfinite( El1 ) ) &
						 np.all( np.isfinite( El2 ) ) &
						 np.all( np.isfinite( El3 ) ) &
						 np.all( np.isfinite( Ev1 ) ) &
						 np.all( np.isfinite( Ev2 ) ) &
						 np.all( np.isfinite( Ev3 ) ) )

	def test_FiniteValues_ini6Cmps(self):
		"""
		Testing if ini6Cmps(shp3) returns FINITE values (no nan/inf allowed)
		"""
		maxshp = 16
		shp3 = np.array(np.random.random_integers(maxshp,high=None,size=3))

		C11,C12,C13,C22,C23,C33 = mpT3libs.ini6Cmps( shp3 )

		self.assertTrue( np.all( np.isfinite( C11 ) ) &
						 np.all( np.isfinite( C12 ) ) &
						 np.all( np.isfinite( C13 ) ) &
						 np.all( np.isfinite( C22 ) ) &
						 np.all( np.isfinite( C23 ) ) &
						 np.all( np.isfinite( C33 ) ) )

	def test_FiniteValues_ini6CmpsLE(self):
		"""
		Testing if ini6CmpsLE(shp3) returns FINITE values (no nan/inf allowed)
		"""
		maxshp = 16
		shp3 = np.array(np.random.random_integers(maxshp,high=None,size=3))

		C11LE,C12LE,C13LE,C22LE,C23LE,C33LE = mpT3libs.ini6CmpsLE( shp3 )

		self.assertTrue( np.all( np.isfinite( C11LE ) ) &
						 np.all( np.isfinite( C12LE ) ) &
						 np.all( np.isfinite( C13LE ) ) &
						 np.all( np.isfinite( C22LE ) ) &
						 np.all( np.isfinite( C23LE ) ) &
						 np.all( np.isfinite( C33LE ) ) )

	def test_FiniteValues_rgl6Cmps(self):
		"""
		Testing if rgl6Cmps(C11,C12,C13,C22,C23,C33) returns FINITE values (no nan/inf allowed)
		Testing if type is np.single (float32)
		"""
		maxshp = 16
		shp3 = np.array(np.random.random_integers(maxshp,high=None,size=3))
		
		C11,C12,C13,C22,C23,C33 = mpT3libs.H3_rnd( shp3 )

		C11,C12,C13,C22,C23,C33 = mpT3libs.rgl6Cmps(C11,C12,C13,C22,C23,C33)

		self.assertTrue( np.all( np.isfinite( C11 ) ) &
						 np.all( np.isfinite( C12 ) ) &
						 np.all( np.isfinite( C13 ) ) &
						 np.all( np.isfinite( C22 ) ) &
						 np.all( np.isfinite( C23 ) ) &
						 np.all( np.isfinite( C33 ) ) )
	
	def test_FiniteValues_mngIdxs_full(self):
		"""
		Testing if mngIdxs(idx,shp3) returns FINITE values (no nan/inf allowed)
		Testing if type idx and numel have negative values (must be both positive)
		"""
		maxshp = 16
		shp3 = np.array(np.random.random_integers(maxshp,high=None,size=3))

		msk3_full = np.ones( shp3 , dtype=bool )

		idx_full = mpT3libs.msk3_to_idx( msk3_full )

		(idx_full,idx_full_ptr,
		 numel_full,numel_full_ptr) = mpT3libs.mngIdxs( idx_full , shp3 )
		
		self.assertTrue( np.all( np.isfinite( idx_full ) ) &
						 np.all( np.isfinite( numel_full ) ) &
						 np.all( np.prod(np.shape(idx_full)) == numel_full ) &
						 np.all( numel_full == np.sum(msk3_full) ) )
		
	def test_FiniteValues_mngIdxs_part(self):
		"""
		Testing if mngIdxs(idx,shp3) returns FINITE values (no nan/inf allowed)
		Testing if type idx and numel have negative values (must be both positive)
		"""
		maxshp = 16
		shp3 = np.array(np.random.random_integers(maxshp,high=None,size=3))

		msk3_part = np.random.rand( shp3[0], shp3[1], shp3[2] ) < 0.5
		
		idx_part = mpT3libs.msk3_to_idx( msk3_part )
		
		(idx_part,idx_part_ptr,
		 numel_part,numel_part_ptr) = mpT3libs.mngIdxs( idx_part , shp3 )
		
		self.assertTrue( np.all( np.isfinite( idx_part ) ) &
						 np.all( np.isfinite( numel_part ) ) &
						 np.all( np.prod(np.shape(idx_part)) == numel_part ) &
						 np.all( numel_part == np.sum(msk3_part) ) )

	def test_FiniteValues_mngIdxs_null(self):
		"""
		Testing if mngIdxs(idx,shp3) returns FINITE values (no nan/inf allowed)
		Testing if type idx and numel have negative values (must be both positive)
		Testing if input idx is empty -> check numel must be 0
		"""
		maxshp = 16
		shp3 = np.array(np.random.random_integers(maxshp,high=None,size=3))

		msk3_null = np.zeros( shp3 , dtype=bool )

		idx_null = mpT3libs.msk3_to_idx( msk3_null )
		
		(idx_null,idx_null_ptr,
		 numel_null,numel_null_ptr) = mpT3libs.mngIdxs( idx_null , shp3 )

		self.assertTrue( np.all( np.isfinite( idx_null ) ) &
						 np.all( np.isfinite( numel_null ) ) &
						 np.all( np.prod(np.shape(idx_null)) == numel_null ) &
						 np.all( numel_null == np.sum(msk3_null) ) )

	def test_FiniteValues_mngIdxs_all(self):
		"""
		Testing if mngIdxs(idx,shp3) returns FINITE values (no nan/inf allowed)
		Testing if type idx and numel have negative values (must be both positive)
		"""
		maxshp = 16
		shp3 = np.array(np.random.random_integers(maxshp,high=None,size=3))

		idx_all = -1

		msk3_all = np.ones( shp3 , dtype=bool )

		(idx_all,idx_all_ptr,
		 numel_all,numel_all_ptr) = mpT3libs.mngIdxs( idx_all , shp3 )

		self.assertTrue( np.all( np.isfinite( idx_all ) ) &
						 np.all( np.isfinite( numel_all ) ) &
						 np.all( np.prod(np.shape(idx_all)) == numel_all ) &
						 np.all( numel_all == np.sum(msk3_all) ) )

	def test_FiniteValues_mpT3_eig(self):
		"""
		Testing if mpT3_eig(H11,H12,H13,H22,H23,H33,idx) returns FINITE values (no nan/inf allowed)
		"""
		maxshp = 16
		shp3 = np.array(np.random.random_integers(maxshp,high=None,size=3))

		H11,H12,H13,H22,H23,H33 = mpT3libs.H3_rnd( shp3 )

		idx_all = -1

		El1,El2,El3,Ev1,Ev2,Ev3,_ = mpT3libs.mpT3_eig( H11,H12,H13,H22,H23,H33, idx_all )

		self.assertTrue( np.all( np.isfinite( El1 ) ) &
						 np.all( np.isfinite( El2 ) ) &
						 np.all( np.isfinite( El3 ) ) &
						 np.all( np.isfinite( Ev1 ) ) &
						 np.all( np.isfinite( Ev2 ) ) &
						 np.all( np.isfinite( Ev3 ) ) )

	def test_FiniteValues_mpT3_to_T3LIC(self):
		"""
		Testing if mpT3_to_T3LIC(El1,El2,El3,Ev1,Ev2,Ev3,idx) returns FINITE values (no nan/inf allowed)
		"""
		maxshp = 16
		shp3 = np.array(np.random.random_integers(maxshp,high=None,size=3))

		El1,El2,El3,Ev1,Ev2,Ev3 = mpT3libs.T3_rnd( shp3 )

		idx_all = -1

		T11,T12,T13,T22,T23,T33 = mpT3libs.mpT3_to_T3LIC( El1,El2,El3,Ev1,Ev2,Ev3, idx_all )

		self.assertTrue( np.all( np.isfinite( T11 ) ) &
						 np.all( np.isfinite( T12 ) ) &
						 np.all( np.isfinite( T13 ) ) &
						 np.all( np.isfinite( T22 ) ) &
						 np.all( np.isfinite( T23 ) ) &
						 np.all( np.isfinite( T33 ) ) )

	def test_FiniteValues_mpT3_to_T3LE(self):
		"""
		Testing if mpT3_to_T3LE(El1,El2,El3,Ev1,Ev2,Ev3,idx) returns FINITE values (no nan/inf allowed)
		"""
		maxshp = 16
		shp3 = np.array(np.random.random_integers(maxshp,high=None,size=3))

		El1,El2,El3,Ev1,Ev2,Ev3 = mpT3libs.T3_rnd( shp3 )

		idx_all = -1

		T11LE,T12LE,T13LE,T22LE,T23LE,T33LE = mpT3libs.mpT3_to_T3LIC( El1,El2,El3,Ev1,Ev2,Ev3, idx_all )

		self.assertTrue( np.all( np.isfinite( T11LE ) ) &
						 np.all( np.isfinite( T12LE ) ) &
						 np.all( np.isfinite( T13LE ) ) &
						 np.all( np.isfinite( T22LE ) ) &
						 np.all( np.isfinite( T23LE ) ) &
						 np.all( np.isfinite( T33LE ) ) )

	def test_FiniteValues_mpT3LIC_to_T3(self):
		"""
		Testing if mpT3LIC_to_T3(T11,T12,T13,T22,T23,T33,idx) returns FINITE values (no nan/inf allowed)
		"""
		maxshp = 16
		shp3 = np.array(np.random.random_integers(maxshp,high=None,size=3))

		T11,T12,T13,T22,T23,T33 = mpT3libs.H3_rnd( shp3 )

		idx_all = -1

		El1,El2,El3,Ev1,Ev2,Ev3,_ = mpT3libs.mpT3LIC_to_T3( T11,T12,T13,T22,T23,T33, idx_all )

		self.assertTrue( np.all( np.isfinite( El1 ) ) &
						 np.all( np.isfinite( El2 ) ) &
						 np.all( np.isfinite( El3 ) ) &
						 np.all( np.isfinite( Ev1 ) ) &
						 np.all( np.isfinite( Ev2 ) ) &
						 np.all( np.isfinite( Ev3 ) ) )

	def test_FiniteValues_mpT3LE_to_T3(self):
		"""
		Testing if mpT3LE_to_T3(T11LE,T12LE,T13LE,T22LE,T23LE,T33LE,idx) returns FINITE values (no nan/inf allowed)
		"""
		maxshp = 16
		shp3 = np.array(np.random.random_integers(maxshp,high=None,size=3))

		T11LE,T12LE,T13LE,T22LE,T23LE,T33LE = mpT3libs.H3_rnd( shp3 )

		idx_all = -1

		El1,El2,El3,Ev1,Ev2,Ev3,_ = mpT3libs.mpT3LE_to_T3( T11LE,T12LE,T13LE,T22LE,T23LE,T33LE, idx_all )

		self.assertTrue( np.all( np.isfinite( El1 ) ) &
						 np.all( np.isfinite( El2 ) ) &
						 np.all( np.isfinite( El3 ) ) &
						 np.all( np.isfinite( Ev1 ) ) &
						 np.all( np.isfinite( Ev2 ) ) &
						 np.all( np.isfinite( Ev3 ) ) )

	# OUTPUT FUNCTIONAL CORRECTNESS TESTS
	def test_ClikeIndexing(self):
		"""
		Testing the C-like Indexing:
			1) Determine a random Hessian or 6 Tensorial Components
			2) Perform Eigen-decomposition for *all* voxels
			3) Perform Eigen-decomposition for *selected* voxels
			4) Compare Eigen-decomposition values:
			4a) Those at the selected voxels MUST match those in the full arrays @ the SAME locations
			4b) Assert that the values for the NON-selected voxels are identically equal to ISOTROPIC
		"""
		maxshp = 16
		shp3 = np.array(np.random.random_integers(maxshp,high=None,size=3))

		errtol = 1e-3

		H11,H12,H13,H22,H23,H33 = mpT3libs.H3_rnd( shp3 )

		msk3_sel = np.random.rand( shp3[0], shp3[1], shp3[2] ) < 0.5
		idx_sel = mpT3libs.msk3_to_idx( msk3_sel )

		idx_all = -1

		El1all,El2all,El3all,Ev1all,Ev2all,Ev3all,_ = mpT3libs.mpT3_eig(H11,H12,H13,H22,H23,H33, idx_all )

		El1sel,El2sel,El3sel,Ev1sel,Ev2sel,Ev3sel,_ = mpT3libs.mpT3_eig(H11,H12,H13,H22,H23,H33, idx_sel )

		El1err = np.array( abs(El1all[msk3_sel] - El1sel[msk3_sel]) )
		El2err = np.array( abs(El2all[msk3_sel] - El2sel[msk3_sel]) )
		El3err = np.array( abs(El3all[msk3_sel] - El3sel[msk3_sel]) )
		Ev1err = np.array( abs(Ev1all[msk3_sel] - Ev1sel[msk3_sel]) )
		Ev2err = np.array( abs(Ev2all[msk3_sel] - Ev2sel[msk3_sel]) )
		Ev3err = np.array( abs(Ev3all[msk3_sel] - Ev3sel[msk3_sel]) )

		El1iso = np.array( El1sel[~msk3_sel] )
		El2iso = np.array( El2sel[~msk3_sel] )
		El3iso = np.array( El3sel[~msk3_sel] )
		
		Ev11iso = np.array( Ev1sel[:,:,:,0] )
		Ev11iso = np.array( Ev11iso[~msk3_sel] )
		Ev12iso = np.array( Ev1sel[:,:,:,1] )
		Ev12iso = np.array( Ev12iso[~msk3_sel] )
		Ev13iso = np.array( Ev1sel[:,:,:,2] )
		Ev13iso = np.array( Ev13iso[~msk3_sel] )

		Ev21iso = np.array( Ev2sel[:,:,:,0] )
		Ev21iso = np.array( Ev21iso[~msk3_sel] )
		Ev22iso = np.array( Ev2sel[:,:,:,1] )
		Ev22iso = np.array( Ev22iso[~msk3_sel] )
		Ev23iso = np.array( Ev2sel[:,:,:,2] )
		Ev23iso = np.array( Ev23iso[~msk3_sel] )

		Ev31iso = np.array( Ev3sel[:,:,:,0] )
		Ev31iso = np.array( Ev31iso[~msk3_sel] )
		Ev32iso = np.array( Ev3sel[:,:,:,1] )
		Ev32iso = np.array( Ev32iso[~msk3_sel] )
		Ev33iso = np.array( Ev3sel[:,:,:,2] )
		Ev33iso = np.array( Ev33iso[~msk3_sel] )

		self.assertTrue( np.all( El1err < errtol )  &
						 np.all( El2err < errtol )  &
						 np.all( El3err < errtol )  &
						 np.all( Ev1err < errtol )  &
						 np.all( Ev2err < errtol )  &
						 np.all( Ev3err < errtol )  &
						 
						 np.all( Ev11iso == 1.0 ) & 
						 np.all( Ev12iso == 0.0 ) &
						 np.all( Ev13iso == 0.0 ) & 

						 np.all( Ev21iso == 0.0 ) & 
						 np.all( Ev22iso == 1.0 ) &
						 np.all( Ev23iso == 0.0 ) & 

						 np.all( Ev31iso == 0.0 ) & 
						 np.all( Ev32iso == 0.0 ) &
						 np.all( Ev33iso == 1.0 ) )

	def test_T3_to_T3LIC_RoundConversion_all(self):
		"""
		Testing the forward-backward Conversion: T3_to_T3LIC -> T3LIC_to_T3
		Input and Output must be equal up to a numeric error.
			1) Determine a random Tensor field
			2a) Perform T3_to_T3LIC conversion for *all* voxels
			2b) Perform T3LIC_to_T3 conversion for *all* voxels
			3) Compare Eigen-decomposition values:
			3a) Compare: OUTPUT of 2b) with INPUT of 2a), for *all* voxels
		"""
		maxshp = 16
		shp3 = np.array(np.random.random_integers(maxshp,high=None,size=3))

		RMSEtol = 1e-3

		El1_in,El2_in,El3_in,Ev1_in,Ev2_in,Ev3_in = mpT3libs.T3_rnd( shp3 )

		idx_all = -1

		T11,T12,T13,T22,T23,T33 = mpT3libs.mpT3_to_T3LIC( El1_in,El2_in,El3_in,Ev1_in,Ev2_in,Ev3_in, idx_all )

		El1_out,El2_out,El3_out,Ev1_out,Ev2_out,Ev3_out,msk3Valid = mpT3libs.mpT3LIC_to_T3( T11,T12,T13,T22,T23,T33, idx_all )

		msk4Valid = np.tile(msk3Valid.reshape( (msk3Valid.shape[0],
												msk3Valid.shape[1],
												msk3Valid.shape[2],
												1                  ) ),3)

		El1err = np.array( abs( El1_in - El1_out ) )
		El2err = np.array( abs( El2_in - El2_out ) )
		El3err = np.array( abs( El3_in - El3_out ) )
		#Ev1err = np.array( abs( abs(Ev1_in) - abs(Ev1_out) ) )
		#Ev2err = np.array( abs( abs(Ev2_in) - abs(Ev2_out) ) )
		#Ev3err = np.array( abs( abs(Ev3_in) - abs(Ev3_out) ) )

		RMSE_El1 = np.power( np.sum( np.power(El1err[msk3Valid],2) ) / np.sum(msk3Valid)  , 0.5)
		RMSE_El2 = np.power( np.sum( np.power(El2err[msk3Valid],2) ) / np.sum(msk3Valid)  , 0.5)
		RMSE_El3 = np.power( np.sum( np.power(El3err[msk3Valid],2) ) / np.sum(msk3Valid)  , 0.5)
		#RMSE_Ev1 = np.power( np.sum( np.power(Ev1err[msk4Valid],2) ) / np.sum(msk4Valid)  , 0.5)
		#RMSE_Ev2 = np.power( np.sum( np.power(Ev2err[msk4Valid],2) ) / np.sum(msk4Valid)  , 0.5)
		#RMSE_Ev3 = np.power( np.sum( np.power(Ev3err[msk4Valid],2) ) / np.sum(msk4Valid)  , 0.5)

		self.assertTrue( np.all( RMSE_El1 < RMSEtol )  &
						 np.all( RMSE_El2 < RMSEtol )  &
						 np.all( RMSE_El3 < RMSEtol )  ) #&
						 #np.all( RMSE_Ev1 < RMSEtol )  &
						 #np.all( RMSE_Ev2 < RMSEtol )  &
						 #np.all( RMSE_Ev3 < RMSEtol )  )
		
	def test_T3_to_T3LIC_RoundConversion_sel(self):
		"""
		Testing the forward-backward Conversion: T3_to_T3LIC -> T3LIC_to_T3
		Input and Output must be equal up to a numeric error.
			1) Determine a random Tensor field
			2a) Perform T3_to_T3LIC conversion for *selected* voxels
			2b) Perform T3LIC_to_T3 conversion for *selected* voxels
			3) Compare Eigen-decomposition values:
			3a) Compare: OUTPUT of 2b) with INPUT of 2a), for *selected* voxels
			3b) Assert:  OUTPUT of 2b) is different from 2a) and is ISOTROPIC, for NON-selected voxels
		"""
		maxshp = 16
		shp3 = np.array(np.random.random_integers(maxshp,high=None,size=3))

		RMSEtol = 1e-3

		El1_in,El2_in,El3_in,Ev1_in,Ev2_in,Ev3_in = mpT3libs.T3_rnd( shp3 )

		msk3_sel = np.random.rand( shp3[0], shp3[1], shp3[2] ) < 0.5
		idx_sel = mpT3libs.msk3_to_idx( msk3_sel )

		T11_sel,T12_sel,T13_sel,T22_sel,T23_sel,T33_sel = mpT3libs.mpT3_to_T3LIC( El1_in,El2_in,El3_in,Ev1_in,Ev2_in,Ev3_in, idx_sel )

		El1_sel,El2_sel,El3_sel,Ev1_sel,Ev2_sel,Ev3_sel,msk3Valid_sel = mpT3libs.mpT3LIC_to_T3( T11_sel,T12_sel,T13_sel,T22_sel,T23_sel,T33_sel, idx_sel )

		msk3CHK = np.logical_and(msk3_sel,msk3Valid_sel)
		#msk4CHK = np.tile( msk3CHK.reshape( (msk3CHK.shape[0],msk3CHK.shape[1],msk3CHK.shape[2], 1 ) ) , 3 )

		El1err = np.array( abs( El1_in - El1_sel ) )
		El2err = np.array( abs( El2_in - El2_sel ) )
		El3err = np.array( abs( El3_in - El3_sel ) )
		#Ev1err = np.array( abs( abs(Ev1_in) - abs(Ev1_sel) ) )
		#Ev2err = np.array( abs( abs(Ev2_in) - abs(Ev2_sel) ) )
		#Ev3err = np.array( abs( abs(Ev3_in) - abs(Ev3_sel) ) )

		RMSE_El1 = np.power( np.sum( np.power(El1err[msk3CHK],2) ) / np.sum(msk3CHK)  , 0.5)
		RMSE_El2 = np.power( np.sum( np.power(El2err[msk3CHK],2) ) / np.sum(msk3CHK)  , 0.5)
		RMSE_El3 = np.power( np.sum( np.power(El3err[msk3CHK],2) ) / np.sum(msk3CHK)  , 0.5)
		#RMSE_Ev1 = np.power( np.sum( np.power(Ev1err[msk4CHK],2) ) / np.sum(msk4CHK)  , 0.5)
		#RMSE_Ev2 = np.power( np.sum( np.power(Ev2err[msk4CHK],2) ) / np.sum(msk4CHK)  , 0.5)
		#RMSE_Ev3 = np.power( np.sum( np.power(Ev3err[msk4CHK],2) ) / np.sum(msk4CHK)  , 0.5)

		self.assertTrue( np.all( RMSE_El1 < RMSEtol )  &
						 np.all( RMSE_El2 < RMSEtol )  &
						 np.all( RMSE_El3 < RMSEtol )  &
						 np.all(El1_sel[~msk3_sel] == 1.0) &
						 np.all(El2_sel[~msk3_sel] == 1.0) &
						 np.all(El3_sel[~msk3_sel] == 1.0) )

		"""
		self.assertTrue( np.all( RMSE_El1 < RMSEtol )  &
						 np.all( RMSE_El2 < RMSEtol )  &
						 np.all( RMSE_El3 < RMSEtol )  &
						 np.all( RMSE_Ev1 < RMSEtol )  &
						 np.all( RMSE_Ev2 < RMSEtol )  &
						 np.all( RMSE_Ev3 < RMSEtol )  &
						 np.all(El1_sel[~msk3_sel] == 1.0) &
						 np.all(El2_sel[~msk3_sel] == 1.0) &
						 np.all(El3_sel[~msk3_sel] == 1.0) )
		"""

	def test_T3LIC_to_T3_RoundConversion_all(self):
		"""
		Testing the forward-backward Conversion: T3LIC_to_T3 -> T3_to_T3LIC
		Input and Output must be equal up to a numeric error.
			1) Determine a random Hessian (or 3D Tensor as 6 Indep. Comps)
			2a) Perform T3LIC_to_T3 conversion for *all* voxels
			2b) Perform T3_to_T3LIC conversion for *all* voxels
			3) Compare Eigen-decomposition values:
			3a) Compare: OUTPUT of 2b) with INPUT of 2a), for *all* voxels
		"""
		maxshp = 16
		shp3 = np.array(np.random.random_integers(maxshp,high=None,size=3))

		RMSEtol = 1e-3

		H11_in,H12_in,H13_in,H22_in,H23_in,H33_in = mpT3libs.H3_rnd( shp3 )
		
		idx_all = -1

		El1,El2,El3,Ev1,Ev2,Ev3,msk3Valid = mpT3libs.mpT3LIC_to_T3( H11_in,H12_in,H13_in,H22_in,H23_in,H33_in, idx_all )

		H11_out,H12_out,H13_out,H22_out,H23_out,H33_out = mpT3libs.mpT3_to_T3LIC( El1,El2,El3,Ev1,Ev2,Ev3, idx_all )


		H11err = np.array( abs( H11_in - H11_out ) )
		H12err = np.array( abs( H12_in - H12_out ) )
		H13err = np.array( abs( H13_in - H13_out ) )
		H22err = np.array( abs( H22_in - H22_out ) )
		H23err = np.array( abs( H23_in - H23_out ) )
		H33err = np.array( abs( H33_in - H33_out ) )

		RMSE_H11 = np.power( np.sum( np.power(H11err[msk3Valid],2) ) / np.sum(msk3Valid)  , 0.5)
		RMSE_H12 = np.power( np.sum( np.power(H12err[msk3Valid],2) ) / np.sum(msk3Valid)  , 0.5)
		RMSE_H13 = np.power( np.sum( np.power(H13err[msk3Valid],2) ) / np.sum(msk3Valid)  , 0.5)
		RMSE_H22 = np.power( np.sum( np.power(H22err[msk3Valid],2) ) / np.sum(msk3Valid)  , 0.5)
		RMSE_H23 = np.power( np.sum( np.power(H23err[msk3Valid],2) ) / np.sum(msk3Valid)  , 0.5)
		RMSE_H33 = np.power( np.sum( np.power(H33err[msk3Valid],2) ) / np.sum(msk3Valid)  , 0.5)

		self.assertTrue( np.all( RMSE_H11 < RMSEtol )  &
						 np.all( RMSE_H12 < RMSEtol )  &
						 np.all( RMSE_H13 < RMSEtol )  &
						 np.all( RMSE_H22 < RMSEtol )  &
						 np.all( RMSE_H23 < RMSEtol )  &
						 np.all( RMSE_H33 < RMSEtol )  )

	def test_T3LIC_to_T3_RoundConversion_sel(self):
		"""
		Testing the forward-backward Conversion: T3LIC_to_T3 -> T3_to_T3LIC
		Input and Output must be equal up to a numeric error.
			1) Determine a random Hessian (or 3D Tensor as 6 Indep. Comps)
			2a) Perform T3_to_T3LIC conversion for *selected* voxels
			2b) Perform T3LIC_to_T3 conversion for *selected* voxels
			3) Compare Eigen-decomposition values:
			3a) Compare: OUTPUT of 2b) with INPUT of 2a), for *selected* voxels
			3b) Assert:  OUTPUT of 2b) is different from 2a) and is ISOTROPIC, for NON-selected voxels
		"""
		maxshp = 16
		shp3 = np.array(np.random.random_integers(maxshp,high=None,size=3))

		RMSEtol = 1e-3

		H11_in,H12_in,H13_in,H22_in,H23_in,H33_in = mpT3libs.H3_rnd( shp3 )
		
		msk3_sel = np.random.rand( shp3[0], shp3[1], shp3[2] ) < 0.5
		idx_sel = mpT3libs.msk3_to_idx( msk3_sel )

		El1,El2,El3,Ev1,Ev2,Ev3,msk3Valid = mpT3libs.mpT3LIC_to_T3( H11_in,H12_in,H13_in,H22_in,H23_in,H33_in, idx_sel )

		H11_out,H12_out,H13_out,H22_out,H23_out,H33_out = mpT3libs.mpT3_to_T3LIC( El1,El2,El3,Ev1,Ev2,Ev3, idx_sel )

		msk3CHK = np.logical_and(msk3_sel,msk3Valid)

		H11err = np.array( abs( H11_in - H11_out ) )
		H12err = np.array( abs( H12_in - H12_out ) )
		H13err = np.array( abs( H13_in - H13_out ) )
		H22err = np.array( abs( H22_in - H22_out ) )
		H23err = np.array( abs( H23_in - H23_out ) )
		H33err = np.array( abs( H33_in - H33_out ) )

		RMSE_H11 = np.power( np.sum( np.power(H11err[msk3CHK],2) ) / np.sum(msk3CHK)  , 0.5)
		RMSE_H12 = np.power( np.sum( np.power(H12err[msk3CHK],2) ) / np.sum(msk3CHK)  , 0.5)
		RMSE_H13 = np.power( np.sum( np.power(H13err[msk3CHK],2) ) / np.sum(msk3CHK)  , 0.5)
		RMSE_H22 = np.power( np.sum( np.power(H22err[msk3CHK],2) ) / np.sum(msk3CHK)  , 0.5)
		RMSE_H23 = np.power( np.sum( np.power(H23err[msk3CHK],2) ) / np.sum(msk3CHK)  , 0.5)
		RMSE_H33 = np.power( np.sum( np.power(H33err[msk3CHK],2) ) / np.sum(msk3CHK)  , 0.5)

		self.assertTrue( np.all( RMSE_H11 < RMSEtol )  &
						 np.all( RMSE_H12 < RMSEtol )  &
						 np.all( RMSE_H13 < RMSEtol )  &
						 np.all( RMSE_H22 < RMSEtol )  &
						 np.all( RMSE_H23 < RMSEtol )  &
						 np.all( RMSE_H33 < RMSEtol )  &
						 np.all( H11_out[~msk3_sel] == 1.0 ) & 
						 np.all( H12_out[~msk3_sel] == 0.0 ) & 
						 np.all( H13_out[~msk3_sel] == 0.0 ) & 
						 np.all( H22_out[~msk3_sel] == 1.0 ) & 
						 np.all( H23_out[~msk3_sel] == 0.0 ) & 
						 np.all( H33_out[~msk3_sel] == 1.0 ) )

	def test_T3_to_T3LE_RoundConversion_all(self):
		"""
		Testing the forward-backward Conversion: T3_to_T3LE -> T3LE_to_T3
		Input and Output must be equal up to a numeric error.
			1) Determine a random 3D Tensor field (Euclidean)
			2a) Perform T3_to_T3LE conversion for *all* voxels
			2b) Perform T3LE_to_T3 conversion for *all* voxels
			3) Compare Eigen-decomposition values:
			3a) Compare: OUTPUT of 2b) with INPUT of 2a), for *all* voxels
		"""
		maxshp = 16
		shp3 = np.array(np.random.random_integers(maxshp,high=None,size=3))

		RMSEtol = 1e-3

		El1_in,El2_in,El3_in,Ev1_in,Ev2_in,Ev3_in = mpT3libs.T3_rnd( shp3 )

		idx_all = -1

		T11LE,T12LE,T13LE,T22LE,T23LE,T33LE = mpT3libs.mpT3_to_T3LE( El1_in,El2_in,El3_in,Ev1_in,Ev2_in,Ev3_in, idx_all )

		El1_out,El2_out,El3_out,Ev1_out,Ev2_out,Ev3_out,msk3Valid = mpT3libs.mpT3LE_to_T3( T11LE,T12LE,T13LE,T22LE,T23LE,T33LE, idx_all )

		#msk4Valid = np.tile(msk3Valid.reshape( (msk3Valid.shape[0],
		#										msk3Valid.shape[1],
		#										msk3Valid.shape[2],
		#										1                  ) ),3)

		El1err = np.array( abs( El1_in - El1_out ) )
		El2err = np.array( abs( El2_in - El2_out ) )
		El3err = np.array( abs( El3_in - El3_out ) )
		#Ev1err = np.array( abs( abs(Ev1_in) - abs(Ev1_out) ) )
		#Ev2err = np.array( abs( abs(Ev2_in) - abs(Ev2_out) ) )
		#Ev3err = np.array( abs( abs(Ev3_in) - abs(Ev3_out) ) )

		RMSE_El1 = np.power( np.sum( np.power(El1err[msk3Valid],2) ) / np.sum(msk3Valid)  , 0.5)
		RMSE_El2 = np.power( np.sum( np.power(El2err[msk3Valid],2) ) / np.sum(msk3Valid)  , 0.5)
		RMSE_El3 = np.power( np.sum( np.power(El3err[msk3Valid],2) ) / np.sum(msk3Valid)  , 0.5)
		#RMSE_Ev1 = np.power( np.sum( np.power(Ev1err[msk4Valid],2) ) / np.sum(msk4Valid)  , 0.5)
		#RMSE_Ev2 = np.power( np.sum( np.power(Ev2err[msk4Valid],2) ) / np.sum(msk4Valid)  , 0.5)
		#RMSE_Ev3 = np.power( np.sum( np.power(Ev3err[msk4Valid],2) ) / np.sum(msk4Valid)  , 0.5)

		self.assertTrue( np.all( RMSE_El1 < RMSEtol )  &
						 np.all( RMSE_El2 < RMSEtol )  &
						 np.all( RMSE_El3 < RMSEtol )  ) #&
						 #np.all( RMSE_Ev1 < RMSEtol )  &
						 #np.all( RMSE_Ev2 < RMSEtol )  &
						 #np.all( RMSE_Ev3 < RMSEtol )  )

	def test_T3_to_T3LE_RoundConversion_sel(self):
		"""
		Testing the forward-backward Conversion: T3_to_T3LE -> T3LE_to_T3
		Input and Output must be equal up to a numeric error.
			1) Determine a random 3D Tensor field (Euclidean)
			2a) Perform T3_to_T3LE conversion for *selected* voxels
			2b) Perform T3LE_to_T3 conversion for *selected* voxels
			3) Compare Eigen-decomposition values:
			3a) Compare: OUTPUT of 2b) with INPUT of 2a), for *selected* voxels
			3b) Assert:  OUTPUT of 2b) is different from 2a) and is ISOTROPIC, for NON-selected voxels
		"""
		maxshp = 16
		shp3 = np.array(np.random.random_integers(maxshp,high=None,size=3))

		RMSEtol = 1e-3

		El1_in,El2_in,El3_in,Ev1_in,Ev2_in,Ev3_in = mpT3libs.T3_rnd( shp3 )

		msk3_sel = np.random.rand( shp3[0], shp3[1], shp3[2] ) < 0.5
		idx_sel = mpT3libs.msk3_to_idx( msk3_sel )

		T11LE_sel,T12LE_sel,T13LE_sel,T22LE_sel,T23LE_sel,T33LE_sel = mpT3libs.mpT3_to_T3LE( El1_in,El2_in,El3_in,Ev1_in,Ev2_in,Ev3_in, idx_sel )

		El1_sel,El2_sel,El3_sel,Ev1_sel,Ev2_sel,Ev3_sel,msk3Valid_sel = mpT3libs.mpT3LE_to_T3( T11LE_sel,T12LE_sel,T13LE_sel,T22LE_sel,T23LE_sel,T33LE_sel, idx_sel )

		msk3CHK = np.logical_and(msk3_sel,msk3Valid_sel)
		#msk4CHK = np.tile( msk3CHK.reshape( (msk3CHK.shape[0],msk3CHK.shape[1],msk3CHK.shape[2], 1 ) ) , 3 )

		El1err = np.array( abs( El1_in - El1_sel ) )
		El2err = np.array( abs( El2_in - El2_sel ) )
		El3err = np.array( abs( El3_in - El3_sel ) )
		#Ev1err = np.array( abs( abs(Ev1_in) - abs(Ev1_sel) ) )
		#Ev2err = np.array( abs( abs(Ev2_in) - abs(Ev2_sel) ) )
		#Ev3err = np.array( abs( abs(Ev3_in) - abs(Ev3_sel) ) )

		RMSE_El1 = np.power( np.sum( np.power(El1err[msk3CHK],2) ) / np.sum(msk3CHK)  , 0.5)
		RMSE_El2 = np.power( np.sum( np.power(El2err[msk3CHK],2) ) / np.sum(msk3CHK)  , 0.5)
		RMSE_El3 = np.power( np.sum( np.power(El3err[msk3CHK],2) ) / np.sum(msk3CHK)  , 0.5)
		#RMSE_Ev1 = np.power( np.sum( np.power(Ev1err[msk4CHK],2) ) / np.sum(msk4CHK)  , 0.5)
		#RMSE_Ev2 = np.power( np.sum( np.power(Ev2err[msk4CHK],2) ) / np.sum(msk4CHK)  , 0.5)
		#RMSE_Ev3 = np.power( np.sum( np.power(Ev3err[msk4CHK],2) ) / np.sum(msk4CHK)  , 0.5)

		self.assertTrue( np.all( RMSE_El1 < RMSEtol )  &
						 np.all( RMSE_El2 < RMSEtol )  &
						 np.all( RMSE_El3 < RMSEtol )  &
						 np.all(El1_sel[~msk3_sel] == 1.0) &
						 np.all(El2_sel[~msk3_sel] == 1.0) &
						 np.all(El3_sel[~msk3_sel] == 1.0) )
		"""
		self.assertTrue( np.all( RMSE_El1 < RMSEtol )  &
						 np.all( RMSE_El2 < RMSEtol )  &
						 np.all( RMSE_El3 < RMSEtol )  &
						 np.all( RMSE_Ev1 < RMSEtol )  &
						 np.all( RMSE_Ev2 < RMSEtol )  &
						 np.all( RMSE_Ev3 < RMSEtol )  &
						 np.all(El1_sel[~msk3_sel] == 1.0) &
						 np.all(El2_sel[~msk3_sel] == 1.0) &
						 np.all(El3_sel[~msk3_sel] == 1.0) )
		"""

	def test_T3LE_to_T3_RoundConversion_all(self):
		"""
		Testing the forward-backward Conversion: T3LE_to_T3 -> T3_to_T3LE
		Input and Output must be equal up to a numeric error.
			1) Determine a random 3D Tensor field (LOG-Euclidean)
			2a) Perform T3LE_to_T3 conversion for *all* voxels
			2b) Perform T3_to_T3LE conversion for *all* voxels
			3) Compare Eigen-decomposition values:
			3a) Compare: OUTPUT of 2b) with INPUT of 2a), for *all* voxels
		"""
		maxshp = 16
		shp3 = np.array(np.random.random_integers(maxshp,high=None,size=3))

		RMSEtol = 1e-3

		El1_ref,El2_ref,El3_ref,Ev1_ref,Ev2_ref,Ev3_ref = mpT3libs.T3_rnd( shp3 )
		idx_all = -1
		T11LE_in,T12LE_in,T13LE_in,T22LE_in,T23LE_in,T33LE_in = mpT3libs.mpT3_to_T3LE( El1_ref,El2_ref,El3_ref,Ev1_ref,Ev2_ref,Ev3_ref, idx_all )

		El1,El2,El3,Ev1,Ev2,Ev3,msk3Valid = mpT3libs.mpT3LE_to_T3( T11LE_in,T12LE_in,T13LE_in,T22LE_in,T23LE_in,T33LE_in, idx_all )

		T11LE_out,T12LE_out,T13LE_out,T22LE_out,T23LE_out,T33LE_out = mpT3libs.mpT3_to_T3LE( El1,El2,El3,Ev1,Ev2,Ev3, idx_all )

		T11LEerr = np.array( abs( T11LE_in - T11LE_out ) )
		T12LEerr = np.array( abs( T12LE_in - T12LE_out ) )
		T13LEerr = np.array( abs( T13LE_in - T13LE_out ) )
		T22LEerr = np.array( abs( T22LE_in - T22LE_out ) )
		T23LEerr = np.array( abs( T23LE_in - T23LE_out ) )
		T33LEerr = np.array( abs( T33LE_in - T33LE_out ) )

		RMSE_T11LE = np.power( np.sum( np.power(T11LEerr[msk3Valid],2) ) / np.sum(msk3Valid)  , 0.5)
		RMSE_T12LE = np.power( np.sum( np.power(T12LEerr[msk3Valid],2) ) / np.sum(msk3Valid)  , 0.5)
		RMSE_T13LE = np.power( np.sum( np.power(T13LEerr[msk3Valid],2) ) / np.sum(msk3Valid)  , 0.5)
		RMSE_T22LE = np.power( np.sum( np.power(T22LEerr[msk3Valid],2) ) / np.sum(msk3Valid)  , 0.5)
		RMSE_T23LE = np.power( np.sum( np.power(T23LEerr[msk3Valid],2) ) / np.sum(msk3Valid)  , 0.5)
		RMSE_T33LE = np.power( np.sum( np.power(T33LEerr[msk3Valid],2) ) / np.sum(msk3Valid)  , 0.5)

		self.assertTrue( np.all( RMSE_T11LE < RMSEtol )  &
						 np.all( RMSE_T12LE < RMSEtol )  &
						 np.all( RMSE_T13LE < RMSEtol )  &
						 np.all( RMSE_T22LE < RMSEtol )  &
						 np.all( RMSE_T23LE < RMSEtol )  &
						 np.all( RMSE_T33LE < RMSEtol )  )

	def test_T3LE_to_T3_RoundConversion_sel(self):
		"""
		Testing the forward-backward Conversion: T3LE_to_T3 -> T3_to_T3LE
		Input and Output must be equal up to a numeric error.
			1) Determine a random 3D Tensor field (LOG-Euclidean)
			2a) Perform T3LE_to_T3 conversion for *selected* voxels
			2b) Perform T3_to_T3LE conversion for *selected* voxels
			3) Compare Eigen-decomposition values:
			3a) Compare: OUTPUT of 2b) with INPUT of 2a), for *selected* voxels
			3b) Assert:  OUTPUT of 2b) is different from 2a) and is ISOTROPIC, for NON-selected voxels
		"""
		maxshp = 16
		shp3 = np.array(np.random.random_integers(maxshp,high=None,size=3))

		RMSEtol = 1e-3

		El1_ref,El2_ref,El3_ref,Ev1_ref,Ev2_ref,Ev3_ref = mpT3libs.T3_rnd( shp3 )
		idx_all = -1
		T11LE_in,T12LE_in,T13LE_in,T22LE_in,T23LE_in,T33LE_in = mpT3libs.mpT3_to_T3LE( El1_ref,El2_ref,El3_ref,Ev1_ref,Ev2_ref,Ev3_ref, idx_all )
		
		msk3_sel = np.random.rand( shp3[0], shp3[1], shp3[2] ) < 0.5
		idx_sel = mpT3libs.msk3_to_idx( msk3_sel )

		El1,El2,El3,Ev1,Ev2,Ev3,msk3Valid = mpT3libs.mpT3LE_to_T3( T11LE_in,T12LE_in,T13LE_in,T22LE_in,T23LE_in,T33LE_in, idx_sel )

		T11LE_out,T12LE_out,T13LE_out,T22LE_out,T23LE_out,T33LE_out = mpT3libs.mpT3_to_T3LE( El1,El2,El3,Ev1,Ev2,Ev3, idx_sel )

		msk3CHK = np.logical_and(msk3_sel,msk3Valid)

		T11LEerr = np.array( abs( T11LE_in - T11LE_out ) )
		T12LEerr = np.array( abs( T12LE_in - T12LE_out ) )
		T13LEerr = np.array( abs( T13LE_in - T13LE_out ) )
		T22LEerr = np.array( abs( T22LE_in - T22LE_out ) )
		T23LEerr = np.array( abs( T23LE_in - T23LE_out ) )
		T33LEerr = np.array( abs( T33LE_in - T33LE_out ) )

		RMSE_T11LE = np.power( np.sum( np.power(T11LEerr[msk3CHK],2) ) / np.sum(msk3CHK)  , 0.5)
		RMSE_T12LE = np.power( np.sum( np.power(T12LEerr[msk3CHK],2) ) / np.sum(msk3CHK)  , 0.5)
		RMSE_T13LE = np.power( np.sum( np.power(T13LEerr[msk3CHK],2) ) / np.sum(msk3CHK)  , 0.5)
		RMSE_T22LE = np.power( np.sum( np.power(T22LEerr[msk3CHK],2) ) / np.sum(msk3CHK)  , 0.5)
		RMSE_T23LE = np.power( np.sum( np.power(T23LEerr[msk3CHK],2) ) / np.sum(msk3CHK)  , 0.5)
		RMSE_T33LE = np.power( np.sum( np.power(T33LEerr[msk3CHK],2) ) / np.sum(msk3CHK)  , 0.5)

		self.assertTrue( np.all( RMSE_T11LE < RMSEtol )  &
						 np.all( RMSE_T12LE < RMSEtol )  &
						 np.all( RMSE_T13LE < RMSEtol )  &
						 np.all( RMSE_T22LE < RMSEtol )  &
						 np.all( RMSE_T23LE < RMSEtol )  &
						 np.all( RMSE_T33LE < RMSEtol )  &
						 np.all( T11LE_out[~msk3_sel] == 0.0 ) & 
						 np.all( T12LE_out[~msk3_sel] == 0.0 ) & 
						 np.all( T13LE_out[~msk3_sel] == 0.0 ) & 
						 np.all( T22LE_out[~msk3_sel] == 0.0 ) & 
						 np.all( T23LE_out[~msk3_sel] == 0.0 ) & 
						 np.all( T33LE_out[~msk3_sel] == 0.0 ) )
