def H3_rnd(shp3):
	# Function to generate random Hessian values for a random image in 3D. 
	#
	# Call: (H11,H12,H13,H22,H23,H33) = H3_rnd( shp3 )
	#
	# *Inputs* 
	# shp3: shape of the 3D image as (dim0, dim1, dim2)
	#
	# *Outputs*
	# H11,H12,H13,H22,H23,H33: 3D arrays of shape 'shp3' with random values.
	# 	They represent the 6 independent partial derivatives of a square and symmetric Hessian matrix H
	# 	evaluated on a random 3D image at voxel (or location) [i,j,k]
	# NB: It is assumed:
	#						[ H11[i,j,k]  H12[i,j,k]  H13[i,j,k] ]
	#        H[i,j,k]   =   [ H12[i,j,k]  H22[i,j,k]  H23[i,j,k] ]
	#						[ H13[i,j,k]  H23[i,j,k]  H33[i,j,k] ]
	#
	# The output values follow a canonical normal distribution (mean=0,std=1).

	# Import Libraries
	import numpy as np
	import warnings

	# Checking Inputs
	if np.prod(np.shape(shp3)) != 3:
		raise Exception('Input: "shp3" should be the shape of a 3D image, e.g. (idx0, idx1, idx2). The value of "shp3" was: {}'.format(shp3))
	if ~ np.all( np.isfinite( shp3 ) , axis=-1 ):
		raise Exception('Input: "shp3" should contain all FINITE elements. The value of "shp3" was: {}'.format(shp3))
	if ~ np.all(shp3 == np.array(shp3,np.int32)):
		shp3 = np.array(shp3,np.int32)
		print(" -- Warning: Input argument (shp3) has non-integer elements - Default: round to closest integers.")
	if ~ np.all( np.greater(shp3 , 0) ):
		raise Exception('Input: "shp3" has non-positive number of elements. The value of "shp3" is: {}'.format(shp3))
	
	# Determining outputs
	H11 = np.random.randn(shp3[0],shp3[1],shp3[2])
	H12 = np.random.randn(shp3[0],shp3[1],shp3[2])
	H13 = np.random.randn(shp3[0],shp3[1],shp3[2])
	H22 = np.random.randn(shp3[0],shp3[1],shp3[2])
	H23 = np.random.randn(shp3[0],shp3[1],shp3[2])
	H33 = np.random.randn(shp3[0],shp3[1],shp3[2])

	(H11,H12,H13,H22,H23,H33) = rgl6Cmps(H11,H12,H13,H22,H23,H33)
	
	return H11,H12,H13,H22,H23,H33

def T3_iso(shp3):
	# Function to generate an isotropic tensor field in 3D
	#
	# Call: (El1,El2,El3,Ev1,Ev2,Ev3) = T3_iso( shp3 )
	#
	# *Inputs* 
	# shp3: shape of the 3D image as (dim0, dim1, dim2)
	#
	# *Outputs*
	# El1,El2,El3: 3D arrays with unitary values of shape shp3.
	# 	They represent the 3 Eigen-values associated to an *isotropic* 
	# 	3D image at voxel (or location) [i,j,k]
	# Ev1,Ev2,Ev3: 4D arrays of shape (shp3,3) with unitary and zero values.
	# 	They represent the 3 canonical Eigen-vectors associated to an *isotropic* 
	# 	3D image at voxel (or location) [i,j,k]

	# Import Libraries
	import numpy as np

	# Checking Inputs
	if np.prod(np.shape(shp3)) != 3:
		raise Exception('Input: "shp3" should be the shape of a 3D image, e.g. (idx0, idx1, idx2). The value of "shp3" was: {}'.format(shp3))
	if ~ np.all( np.isfinite( shp3 ) , axis=-1 ):
		raise Exception('Input: "shp3" should contain all FINITE elements. The value of "shp3" was: {}'.format(shp3))
	if ~ np.all(shp3 == np.array(shp3,np.int32)):
		shp3 = np.array(shp3,np.int32)
		print(" -- Warning: Input argument (shp3) has non-integer elements - Default: round to closest integers.")
	if ~ np.all( np.greater(shp3 , 0) ):
		raise Exception('Input: "shp3" has non-positive number of elements. The value of "shp3" is: {}'.format(shp3))

	# Initialising outputs
	(El1,El2,El3,
	 Ev11,Ev12,Ev13,
	 Ev21,Ev22,Ev23,
	 Ev31,Ev32,Ev33) = iniElvs(shp3)

	Ev1,Ev2,Ev3 = catEvs(Ev11,Ev12,Ev13,
						 Ev21,Ev22,Ev23,
						 Ev31,Ev32,Ev33)

	return El1,El2,El3,Ev1,Ev2,Ev3

def T3_rnd(shp3):
	# Function to generate a random tensor field in 3D
	#
	# Call: (El1,El2,El3,Ev1,Ev2,Ev3) = T3_rnd( shp3 )
	#
	# *Inputs*
	# shp3: shape of the 3D image as (dim0, dim1, dim2)
	#
	# *Outputs*
	# El1,El2,El3: 3D arrays with random values of shape shp3,
	# 	they represent the 3 Eigen-values associated to a
	# 	3D image at voxel (or location) [idx0,idx1,idx2]
	# Ev1,Ev2,Ev3: 4D arrays with values of shape (shp3,3),
	# 	they represent the 3 Eigen-vectors associated to the 
	# 	eigen-values of a 3D image at voxel (or location) [idx0,idx1,idx2]
	#
	# *** NB: Eigen-values and associated Eigen-vectors are *SORTED* to satisfy
	#         the condition: abs(El1) <= abs(El2) <= abs(El3)

	# Import Libraries
	import numpy as np

	# Checking Inputs
	if np.prod(np.shape(shp3)) != 3:
		raise Exception('Input: "shp3" should be the shape of a 3D image, e.g. (idx0, idx1, idx2). The value of "shp3" was: {}'.format(shp3))
	if ~ np.all( np.isfinite( shp3 ) , axis=-1 ):
		raise Exception('Input: "shp3" should contain all FINITE elements. The value of "shp3" was: {}'.format(shp3))
	if ~ np.all(shp3 == np.array(shp3,np.int32)):
		shp3 = np.array(shp3,np.int32)
		print(" -- Warning: Input argument (shp3) has non-integer elements - Default: round to closest integers.")
	if ~ np.all( np.greater(shp3 , 0) ):
		raise Exception('Input: "shp3" has non-positive number of elements. The value of "shp3" is: {}'.format(shp3))

	# Generating random Eigen-values
	psi = np.random.rand(shp3[0],shp3[1],shp3[2],3) # psi is 4D!
	psi_divisor = np.power( np.prod(psi,axis=-1) , 1/3.0 )
	psi_divisor = psi_divisor.reshape( ( psi_divisor.shape[0],psi_divisor.shape[1],psi_divisor.shape[2], 1 ) )
	psi = np.divide( psi , np.tile( psi_divisor , 3 ) )

	# Correcting for possible non-finite values (nan/inf)
	psiFINITE = np.all( np.isfinite( psi ) , axis=-1 )
	psiFINITE = psiFINITE.reshape( ( psiFINITE.shape[0] , psiFINITE.shape[1] , psiFINITE.shape[2] , 1) )
	psiFINITE = np.tile( psiFINITE , 3 )
	psi[~psiFINITE] = 1.0

	# Sorting in Ascending Order
	psi = np.sort(psi, axis=-1)

	# Generating 3 random Eigen-vectors
	# First Eigen-vector
	q1 = np.random.rand(shp3[0],shp3[1],shp3[2],3)
	q1_norm = np.linalg.norm( q1 , 2 , axis=-1 )
	q1_norm = q1_norm.reshape( (q1_norm.shape[0],q1_norm.shape[1],q1_norm.shape[2], 1 ) )
	q1 = np.divide( q1 , np.tile( q1_norm , 3 ) )# unit-vector
	# Correcting for possible non-finite values (nan/inf)
	q1FINITE = np.all( np.isfinite( q1 ) , axis=-1 ) # 3D!
	q10 = q1[:,:,:,0] # pointers!
	q11 = q1[:,:,:,1] # pointers!
	q12 = q1[:,:,:,2] # pointers!
	q10[~q1FINITE] = 1.0
	q11[~q1FINITE] = 0.0
	q12[~q1FINITE] = 0.0

	# Second Eigen-vector
	q2 = np.random.rand(shp3[0],shp3[1],shp3[2],3)
	dot_q2q1 = np.sum( np.multiply( q2 , q1 ) , axis=-1 )
	dot_q2q1 = dot_q2q1.reshape( (dot_q2q1.shape[0],dot_q2q1.shape[1],dot_q2q1.shape[2], 1 ) )
	q2 = q2 - ( np.multiply( np.tile( dot_q2q1 , 3 ) , q1 ) )
	q2_norm = np.linalg.norm( q2 , 2 , axis=-1 )
	q2_norm = q2_norm.reshape( (q2_norm.shape[0],q2_norm.shape[1],q2_norm.shape[2], 1 ) )
	q2 = np.divide( q2 , np.tile( q2_norm , 3 ) )# unit-vector
	# Correcting for possible non-finite values (nan/inf)
	q2FINITE = np.all( np.isfinite( q2 ) , axis=-1 ) # 3D!
	q20 = q2[:,:,:,0] # pointers!
	q21 = q2[:,:,:,1] # pointers!
	q22 = q2[:,:,:,2] # pointers!
	q20[~q2FINITE] = 0.0 
	q21[~q2FINITE] = 1.0
	q22[~q2FINITE] = 0.0

	# Third Eigen-vector
	q3 = np.cross( q1 , q2 , axis=-1)
	q3_norm = np.linalg.norm( q3 , 2 , axis=-1 )
	q3_norm = q3_norm.reshape( (q3_norm.shape[0],q3_norm.shape[1],q3_norm.shape[2], 1 ) )
	q3 = np.divide( q3 , np.tile( q3_norm , 3 ) )# unit-vector
	# Correcting for possible non-finite values (nan/inf)
	q3FINITE = np.all( np.isfinite( q3 ) , axis=-1 ) # 3D!
	q30 = q3[:,:,:,0] # pointers!
	q31 = q3[:,:,:,1] # pointers!
	q32 = q3[:,:,:,2] # pointers!
	q30[~q3FINITE] = 0.0 
	q31[~q3FINITE] = 0.0
	q32[~q3FINITE] = 1.0

	El1,El2,El3,Ev1,Ev2,Ev3 = rglElvs( psi[:,:,:,0], psi[:,:,:,1], psi[:,:,:,2], q1, q2, q3 )

	return El1,El2,El3,Ev1,Ev2,Ev3

def msk3_to_idx(msk3):
	# Function to determine the C-like linear indices, given a logical mask of a 3D image.
	#
	# Call: idx = msk3_to_idx( msk3 )
	#
	# *Inputs*
	# msk3: boolean mask of a 3D image of shape (dim0, dim1, dim2)
	#
	# *Outputs*
	# idx: C-like linear indices corresponding to the 'True' voxels in msk3
	#	e.g. if msk3 has shape (3,3,3) and the only true value is in the middle,
	#	i.e. msk3[1,1,1] = True, the idx array will be a scalar equal to 13.

	# Import Libraries
	import numpy as np

	# Checking Inputs
	if np.prod(np.shape(np.shape(msk3))) != 3:
		raise Exception('Input: "msk3" should have the shape of a 3D image, e.g. (idx0, idx1, idx2). The shape was: {}'.format( np.shape(msk3) ))
	if msk3.dtype != bool:
		raise Exception('Input: "msk3" should be logical, i.e. boolean type. The data type was: {}'.format( msk3.dtype ))

	idx = np.array( np.ravel_multi_index( np.nonzero(msk3) , msk3.shape ) , dtype=np.int32)

	return idx

def inimsk3Valid(shp3):
	# Function to initialise the logical validity mask variables used in the compiled Shared Library
	# This function set the data type to boolean, and returns the associated pointer.
	#
	# Call: (msk3Valid,msk3Valid_ptr) = inimsk3Valid( shp3 )
	#
	# *Inputs*
	# shp3: shape of the 3D image as (dim0, dim1, dim2)
	#
	# *Outputs*
	# msk3Valid: logical validity mask of 3D shape (shp3), initialised as 'True'
	# msk3Valid_ptr: pointer of the logical validity mask.

	# Import Libraries
	import numpy as np
	import ctypes

	if np.prod(np.shape(shp3)) != 3:
		raise Exception('Input: "shp3" should be the shape of a 3D image, e.g. (idx0, idx1, idx2). The value of "shp3" was: {}'.format(shp3))
	if ~ np.all( np.isfinite( shp3 ) , axis=-1 ):
		raise Exception('Input: "shp3" should contain all FINITE elements. The value of "shp3" was: {}'.format(shp3))
	if ~ np.all(shp3 == np.array(shp3,np.int32)):
		shp3 = np.array(shp3,np.int32)
		print(" -- Warning: Input argument (shp3) has non-integer elements - Default: round to closest integers.")
	if ~ np.all( np.greater(shp3 , 0) ):
		raise Exception('Input: "shp3" has non-positive number of elements. The value of "shp3" is: {}'.format(shp3))

	msk3Valid = np.ones( shp3 ,dtype=np.bool_)
	msk3Valid_ptr = msk3Valid.ctypes.data_as( ctypes.POINTER( ctypes.c_bool ) )

	return msk3Valid,msk3Valid_ptr

def iniElvs(shp3):
	# Function to initialise the Eigen-Decomposition variables used in the compiled Shared Library
	# This function set the data type to float32 (i.e. np.single).
	#
	# Call: (El1,El2,El3,
	#		 Ev11,Ev12,Ev13,
	#		 Ev21,Ev22,Ev23,
	#		 Ev31,Ev32,Ev33) = iniElvs( shp3 )
	#
	# *Inputs*
	# shp3: shape of the 3D image as (dim0, dim1, dim2)
	#
	# *Outputs*
	# El1,El2,El3: Eigen-values of 3D shape (shp3)
	# Ev11,Ev12,Ev13,...,Ev32,Ev33: Separate components of the 4D Eigen-vectors.
	#	Each output component has the same 3D shape (i.e. shp3).

	# Import Libraries
	import numpy as np

	if np.prod(np.shape(shp3)) != 3:
		raise Exception('Input: "shp3" should be the shape of a 3D image, e.g. (idx0, idx1, idx2). The value of "shp3" was: {}'.format(shp3))
	if ~ np.all( np.isfinite( shp3 ) , axis=-1 ):
		raise Exception('Input: "shp3" should contain all FINITE elements. The value of "shp3" was: {}'.format(shp3))
	if ~ np.all(shp3 == np.array(shp3,np.int32)):
		shp3 = np.array(shp3,np.int32)
		print(" -- Warning: Input argument (shp3) has non-integer elements - Default: round to closest integers.")
	if ~ np.all( np.greater(shp3 , 0) ):
		raise Exception('Input: "shp3" has non-positive number of elements. The value of "shp3" is: {}'.format(shp3))

	# Initialising Outputs
	El1 = np.ones( shp3 ,dtype=np.single)
	El2 = np.ones( shp3 ,dtype=np.single)
	El3 = np.ones( shp3 ,dtype=np.single)
	
	Ev11 = np.ones( shp3 ,dtype=np.single)
	Ev12 = np.zeros( shp3 ,dtype=np.single)
	Ev13 = np.zeros( shp3 ,dtype=np.single)

	Ev21 = np.zeros( shp3 ,dtype=np.single)
	Ev22 = np.ones( shp3 ,dtype=np.single)
	Ev23 = np.zeros( shp3 ,dtype=np.single)
	
	Ev31 = np.zeros( shp3 ,dtype=np.single)
	Ev32 = np.zeros( shp3 ,dtype=np.single)
	Ev33 = np.ones( shp3 ,dtype=np.single)

	return El1,El2,El3,Ev11,Ev12,Ev13,Ev21,Ev22,Ev23,Ev31,Ev32,Ev33

def catEvs(Ev11,Ev12,Ev13,Ev21,Ev22,Ev23,Ev31,Ev32,Ev33):
	# Function to concatenate the Eigen-Vector Components as 4D variables
	#
	# Call: (Ev1,Ev2,Ev3) = catEvs( Ev11,Ev12,Ev13,Ev21,Ev22,Ev23,Ev31,Ev32,Ev33 )
	#
	# *Inputs*
	# Ev11,Ev12,Ev13,...,Ev32,Ev33: Separate components of the 4D Eigen-vectors.
	#	Each input component must have the same 3D shape (i.e. shp3). 
	#
	# *Outputs*
	# Ev1,Ev2,Ev3: Concatenated Eigen-vectors array of 4D shape (shp3,3)

	# Import Libraries
	import numpy as np

	if  ( (np.shape(Ev11) != np.shape(Ev12)) & (np.shape(Ev12) != np.shape(Ev13)) &
		  (np.shape(Ev13) != np.shape(Ev21)) & (np.shape(Ev21) != np.shape(Ev22)) &
		  (np.shape(Ev22) != np.shape(Ev23)) & (np.shape(Ev23) != np.shape(Ev31)) &
		  (np.shape(Ev31) != np.shape(Ev32)) & (np.shape(Ev32) != np.shape(Ev33)) ):
		raise Exception('All inputs should have the *SAME* shape. The shape of the first argument is: {}'.format( np.shape(Ev11) ) )

	shp3 = np.shape(Ev11)

	Ev1 = np.concatenate( ( Ev11.reshape( (shp3[0], shp3[1], shp3[2], 1) ),
						    Ev12.reshape( (shp3[0], shp3[1], shp3[2], 1) ),
						    Ev13.reshape( (shp3[0], shp3[1], shp3[2], 1) ) ), axis=-1)

	Ev2 = np.concatenate( ( Ev21.reshape( (shp3[0], shp3[1], shp3[2], 1) ),
						    Ev22.reshape( (shp3[0], shp3[1], shp3[2], 1) ),
						    Ev23.reshape( (shp3[0], shp3[1], shp3[2], 1) ) ), axis=-1)

	Ev3 = np.concatenate( ( Ev31.reshape( (shp3[0], shp3[1], shp3[2], 1) ),
							Ev32.reshape( (shp3[0], shp3[1], shp3[2], 1) ),
							Ev33.reshape( (shp3[0], shp3[1], shp3[2], 1) ) ), axis=-1)

	return Ev1,Ev2,Ev3

def sepEvs(Ev1,Ev2,Ev3):
	# Function to separate the Eigen-Vector Components from 4D variables to 3D variables
	# This function set the data type to float32 (i.e. np.single).
	#
	# Call: (Ev11,Ev12,Ev13,Ev21,Ev22,Ev23,Ev31,Ev32,Ev33) = sepEvs( Ev1,Ev2,Ev3 )
	#
	# *Inputs*
	# Ev1,Ev2,Ev3: Concatenated Eigen-vectors array of 4D shape (shp3,3).
	#
	# *Outputs*
	# Ev11,Ev12,Ev13,...,Ev32,Ev33: Separate components of the 4D Eigen-vectors.
	#	Each output components has the same 3D shape (i.e. shp3). 

	# Import Libraries
	import numpy as np

	if ( (np.shape(Ev1) != np.shape(Ev2)) & (np.shape(Ev2) != np.shape(Ev3)) ):
		raise Exception('All inputs should have the *SAME* shape. The shape of the first argument is: {}'.format( np.shape(Ev1) ) )
	
	shp4 = np.shape(Ev1)
	
	if ( np.prod(np.shape(shp4)) != 4 ):
		raise Exception('All Eigen-vectors must have the shape of a 4D image. The shape of the first argument is: {}'.format( np.shape(Ev1) ) )
	if ( shp4[-1] != 3 ):
		raise Exception('All Eigen-vectors must have the 3 elements in the last dimension. The number of elements in the last dimension is: {}'.format( shp4[-1] ) )

	Ev11 = np.array( Ev1[:,:,:,0], dtype=np.single);
	Ev12 = np.array( Ev1[:,:,:,1], dtype=np.single);
	Ev13 = np.array( Ev1[:,:,:,2], dtype=np.single);
	
	Ev21 = np.array( Ev2[:,:,:,0], dtype=np.single);
	Ev22 = np.array( Ev2[:,:,:,1], dtype=np.single);
	Ev23 = np.array( Ev2[:,:,:,2], dtype=np.single);

	Ev31 = np.array( Ev3[:,:,:,0], dtype=np.single);
	Ev32 = np.array( Ev3[:,:,:,1], dtype=np.single);
	Ev33 = np.array( Ev3[:,:,:,2], dtype=np.single);

	return Ev11,Ev12,Ev13,Ev21,Ev22,Ev23,Ev31,Ev32,Ev33

def rglElvs(El1,El2,El3,Ev1,Ev2,Ev3):
	# Function to regularise the Eigen-decomposition variables.
	# This function set the data type to float32 (i.e. np.single).
	#
	# Call: (El1,El2,El3,Ev1,Ev2,Ev3) = rglElvs( El1,El2,El3,Ev1,Ev2,Ev3 )
	#
	# *Inputs*
	# El1,El2,El3: Eigen-values of 3D shape (shp3).
	# Ev1,Ev2,Ev3: Eigen-vectors of 4D shape (shp3,3).
	#
	# *Outputs*
	# Same as Inputs, but in float32 data type

	# Import Libraries
	import numpy as np

	if ( (np.shape(El1) != np.shape(El2)) & (np.shape(El2) != np.shape(El3)) ):
		raise Exception('All Eigen-values should have the *SAME* shape. The shape of the first Eigen-value is: {}'.format( np.shape(El1) ) )

	shp3 = np.shape(El1)

	if ( np.prod(np.shape(shp3)) != 3 ):
		raise Exception('All Eigen-values must have the shape of a 3D image. The shape of the first Eigen-value is: {}'.format( np.shape(El1) ) )
	 
	if ( (np.shape(Ev1) != np.shape(Ev2)) & (np.shape(Ev2) != np.shape(Ev3)) ):
		raise Exception('All Eigen-vectors should have the *SAME* shape. The shape of the first Eigen-vector is: {}'.format( np.shape(Ev1) ) )
	
	shp4 = np.shape(Ev1)
	
	if ( np.prod(np.shape(shp4)) != 4 ):
		raise Exception('All Eigen-vectors must have the shape of a 4D image. The shape of the first argument is: {}'.format( np.shape(Ev1) ) )
	
	if ( (shp3[0] != shp4[0]) & (shp3[1] != shp4[1]) & (shp3[2] != shp4[2])):
		raise Exception('Eigen-values and Eigen-vectors shapes do not match the shape of the underlying 3D image. The shapes of the first Eigen-value and Eigen-vector is: {}'.format( [np.shape(El1),np.shape(Ev1)] ) )

	if ( shp4[-1] != 3 ):
		raise Exception('All Eigen-vectors must have 3 elements in the last dimension. The number of elements in the last dimension is: {}'.format( shp4[-1] ) )


	El1_out = np.array( El1 , dtype=np.single )
	El2_out = np.array( El2 , dtype=np.single )
	El3_out = np.array( El3 , dtype=np.single )
	Ev1_out = np.array( Ev1 , dtype=np.single )
	Ev2_out = np.array( Ev2 , dtype=np.single )
	Ev3_out = np.array( Ev3 , dtype=np.single )

	return El1_out,El2_out,El3_out,Ev1_out,Ev2_out,Ev3_out

def ptrElvs(El1,El2,El3,Ev11,Ev12,Ev13,Ev21,Ev22,Ev23,Ev31,Ev32,Ev33):
	# Function to determine the pointers to the Eigen-Decomposition variables used in the compiled Shared Library
	#
	# Call: (El1_ptr,El2_ptr,El3_ptr,
	#   	 Ev11_ptr,Ev12_ptr,Ev13_ptr,
	#		 Ev21_ptr,Ev22_ptr,Ev23_ptr,
	#		 Ev31_ptr,Ev32_ptr,Ev33_ptr) = ptrElvs( El1,El2,El3,
	#											    Ev11,Ev12,Ev13,
	#											    Ev21,Ev22,Ev23,
	#											    Ev31,Ev32,Ev33 )
	#
	# *Inputs*
	# El1,El2,El3: Eigen-values of shape (shp3)
	# Ev11,Ev12, Ev13,...,Ev32,Ev33: Separate components of the 4D Eigen-vectors.
	#	Each input component has the same 3D shape (i.e. shp3). 
	#
	# *Outputs*
	# El1_ptr,El2_ptr,El3_ptr,..., Ev33_ptr: C-like pointers to the input 3D arrays.

	# Import Libraries
	import numpy as np
	import ctypes

	if  ( (np.shape(El1)  != np.shape(El2))  & (np.shape(El2)  != np.shape(El3))  & (np.shape(El3) != np.shape(Ev11)) &
		  (np.shape(Ev11) != np.shape(Ev12)) & (np.shape(Ev12) != np.shape(Ev13)) &
		  (np.shape(Ev13) != np.shape(Ev21)) & (np.shape(Ev21) != np.shape(Ev22)) &
		  (np.shape(Ev22) != np.shape(Ev23)) & (np.shape(Ev23) != np.shape(Ev31)) &
		  (np.shape(Ev31) != np.shape(Ev32)) & (np.shape(Ev32) != np.shape(Ev33)) ):
		raise Exception('All inputs should have the *SAME* shape. The shape of the first argument is: {}'.format( np.shape(El1) ) )

	if ( (El1.dtype  != 'float32') & (El2.dtype  != 'float32') & (El3.dtype  != 'float32') &
		 (Ev11.dtype != 'float32') & (Ev12.dtype != 'float32') & (Ev13.dtype != 'float32') &
		 (Ev21.dtype != 'float32') & (Ev22.dtype != 'float32') & (Ev23.dtype != 'float32') &
		 (Ev31.dtype != 'float32') & (Ev32.dtype != 'float32') & (Ev33.dtype != 'float32') ):
		raise Exception('All inputs must have the same *DATA TYPE* as float32 (i.e. numpy.single). The data type of the first argument is: {}'.format( El1.dtype ) )

	El1_ptr = El1.ctypes.data_as( ctypes.POINTER( ctypes.c_float ) )
	El2_ptr = El2.ctypes.data_as( ctypes.POINTER( ctypes.c_float ) )
	El3_ptr = El3.ctypes.data_as( ctypes.POINTER( ctypes.c_float ) )
	
	Ev11_ptr = Ev11.ctypes.data_as( ctypes.POINTER( ctypes.c_float ) )
	Ev12_ptr = Ev12.ctypes.data_as( ctypes.POINTER( ctypes.c_float ) )
	Ev13_ptr = Ev13.ctypes.data_as( ctypes.POINTER( ctypes.c_float ) )
	
	Ev21_ptr = Ev21.ctypes.data_as( ctypes.POINTER( ctypes.c_float ) )
	Ev22_ptr = Ev22.ctypes.data_as( ctypes.POINTER( ctypes.c_float ) )
	Ev23_ptr = Ev23.ctypes.data_as( ctypes.POINTER( ctypes.c_float ) )
	
	Ev31_ptr = Ev31.ctypes.data_as( ctypes.POINTER( ctypes.c_float ) )
	Ev32_ptr = Ev32.ctypes.data_as( ctypes.POINTER( ctypes.c_float ) )
	Ev33_ptr = Ev33.ctypes.data_as( ctypes.POINTER( ctypes.c_float ) )

	return El1_ptr,El2_ptr,El3_ptr,Ev11_ptr,Ev12_ptr,Ev13_ptr,Ev21_ptr,Ev22_ptr,Ev23_ptr,Ev31_ptr,Ev32_ptr,Ev33_ptr

def ini6Cmps(shp3):
	# Function to initialise the 6 Indep. Component variables used in the compiled Shared Library
	#
	# Call: (C11,C12,C13,C22,C23,C33) = ini6Cmps( shp3 )
	#
	# *Inputs*
	# shp3: shape of the 3D array
	#
	# *Outputs*
	# C11,C12,C13,C22,C23,C33: 3D arrays with initialised values of shape shp3
	# 	they represent the 6 independent components (e.g. partial derivatives of a Hessian matrix, or
	# 	the 6 linearly independent tensorial components in the EUCLIDEAN space, 
	#   evaluated on a 3D image at voxel (or location) [i,j,k]
	# NB: The initialisation sets the outputs as *Isotropic*

	# Import Libraries
	import numpy as np

	# Checking Inputs
	if np.prod(np.shape(shp3)) != 3:
		raise Exception('Input: "shp3" should be the shape of a 3D image, e.g. (idx0, idx1, idx2). The value of "shp3" was: {}'.format(shp3))
	if ~ np.all( np.isfinite( shp3 ) , axis=-1 ):
		raise Exception('Input: "shp3" should contain all FINITE elements. The value of "shp3" was: {}'.format(shp3))
	if ~ np.all(shp3 == np.array(shp3,np.int32)):
		shp3 = np.array(shp3,np.int32)
		print(" -- Warning: Input argument (shp3) has non-integer elements - Default: round to closest integers.")
	if ~ np.all( np.greater(shp3 , 0) ):
		raise Exception('Input: "shp3" has non-positive number of elements. The value of "shp3" is: {}'.format(shp3))

	C11 = np.ones(shp3,dtype=np.single)
	C12 = np.zeros(shp3,dtype=np.single)
	C13 = np.zeros(shp3,dtype=np.single)
	C22 = np.ones(shp3,dtype=np.single)
	C23 = np.zeros(shp3,dtype=np.single)
	C33 = np.ones(shp3,dtype=np.single)

	return C11,C12,C13,C22,C23,C33

def ini6CmpsLE(shp3):
	# Function to initialise the 6 Indep. Component variables in the LOG-EUCLIDEAN space
	# used in the compiled Shared Library
	#
	# Call: (C11LE,C12LE,C13LE,C22LE,C23LE,C33LE) = ini6CmpsLE( shp3 )
	#
	# *Inputs*
	# shp3: shape of the 3D image as (dim0, dim1, dim2)
	#
	# *Outputs*
	# C11LE,C12LE,C13LE,C22LE,C23LE,C33LE: 3D arrays with initialised values of shape shp3
	# 	they represent the 6 linearly independent components of a tensor field in 3D
	#	in the *LOG-EUCLIDEAN* space.
	# NB: The initialisation sets the outputs as *Isotropic* in the LOG-EUCLIDEAN space

	# Import Libraries
	import numpy as np

	# Checking Inputs
	if np.prod(np.shape(shp3)) != 3:
		raise Exception('Input: "shp3" should be the shape of a 3D image, e.g. (idx0, idx1, idx2). The value of "shp3" was: {}'.format(shp3))
	if ~ np.all( np.isfinite( shp3 ) , axis=-1 ):
		raise Exception('Input: "shp3" should contain all FINITE elements. The value of "shp3" was: {}'.format(shp3))
	if ~ np.all(shp3 == np.array(shp3,np.int32)):
		shp3 = np.array(shp3,np.int32)
		print(" -- Warning: Input argument (shp3) has non-integer elements - Default: round to closest integers.")
	if ~ np.all( np.greater(shp3 , 0) ):
		raise Exception('Input: "shp3" has non-positive number of elements. The value of "shp3" is: {}'.format(shp3))

	C11LE = np.zeros(shp3,dtype=np.single)
	C12LE = np.zeros(shp3,dtype=np.single)
	C13LE = np.zeros(shp3,dtype=np.single)
	C22LE = np.zeros(shp3,dtype=np.single)
	C23LE = np.zeros(shp3,dtype=np.single)
	C33LE = np.zeros(shp3,dtype=np.single)

	return C11LE,C12LE,C13LE,C22LE,C23LE,C33LE

def rgl6Cmps(C11,C12,C13,C22,C23,C33):
	# Function to regularise the 6 Indep. Component variables used in the compiled Shared Library.
	# This function set the data type to float32 (i.e. np.single).
	#
	# Call: (C11,C12,C13,C22,C23,C33) = rgl6Cmps( C11,C12,C13,C22,C23,C33 )
	#
	# *Inputs*
	# C11,C12,C13,C22,C23,C33: 3D arrays with 6 independent (C)omponents.
	# 	All C11,C12,...,C33 must have the same shape: shp3 = C11.shape, 
	# 	corresponding to a 3D image for each voxel (or location) at [i,j,k].
	#
	# *Outputs*
	# Same as Inputs, but in float32 data type

	# Import Libraries
	import numpy as np

	if  ( (np.shape(C11) != np.shape(C12)) &
		  (np.shape(C12) != np.shape(C13)) &
		  (np.shape(C13) != np.shape(C22)) &
		  (np.shape(C22) != np.shape(C23)) &
		  (np.shape(C23) != np.shape(C33)) ):
		raise Exception('All inputs should have the *SAME* shape. The shape of the first argument is: {}'.format( np.shape(C11) ) )

	shp3 = np.shape(C11)

	if np.prod(np.shape(shp3)) != 3:
		raise Exception('Inputs should all have the shape of a 3D image, e.g. (idx0, idx1, idx2). The shape was: {}'.format(shp3))

	C11_out = np.array(C11,dtype=np.single)
	C12_out = np.array(C12,dtype=np.single)
	C13_out = np.array(C13,dtype=np.single)
	C22_out = np.array(C22,dtype=np.single)
	C23_out = np.array(C23,dtype=np.single)
	C33_out = np.array(C33,dtype=np.single)

	return C11_out,C12_out,C13_out,C22_out,C23_out,C33_out

def ptr6Cmps(C11,C12,C13,C22,C23,C33):
	# Function to determine the pointers to the 6 Indep. Component variables used in the compiled Shared Library
	#
	# (C11_ptr,C12_ptr,C13_ptr,
	#  C22_ptr,C23_ptr,C33_ptr) = ptr6Cmps( C11,C12,C13,C22,C23,C33 )
	#
	# *Inputs*
	# C11,C12,C13,C22,C23,C33: 3D arrays with 6 independent (C)omponents.
	# 	All C11,C12,...,C33 must have the same shape: shp3 = C11.shape, 
	# 	corresponding to a 3D image for each voxel (or location) at [i,j,k].
	#
	# *Outputs*
	# C11_ptr,C12_ptr,...,C33_ptr: C-like pointers to the 3D arrays with 6 independent components.

	# Import Libraries
	import numpy as np
	import ctypes

	if  ( (np.shape(C11) != np.shape(C12)) &
		  (np.shape(C12) != np.shape(C13)) & 
		  (np.shape(C13) != np.shape(C22)) &
		  (np.shape(C22) != np.shape(C23)) &
		  (np.shape(C23) != np.shape(C33)) ):
		raise Exception('All inputs should have the *SAME* shape. The shape of the first argument is: {}'.format( np.shape(C11) ) )

	if ( (C11.dtype != 'float32') & (C12.dtype != 'float32') & (C13.dtype != 'float32') &
		 (C22.dtype != 'float32') & (C23.dtype != 'float32') & (C33.dtype != 'float32') ):
		raise Exception('All inputs must have the same *DATA TYPE* as float32 (i.e. numpy.single). The data type of the first argument is: {}'.format( C11.dtype ) )

	C11_ptr = C11.ctypes.data_as( ctypes.POINTER( ctypes.c_float ) )
	C12_ptr = C12.ctypes.data_as( ctypes.POINTER( ctypes.c_float ) )
	C13_ptr = C13.ctypes.data_as( ctypes.POINTER( ctypes.c_float ) )
	C22_ptr = C22.ctypes.data_as( ctypes.POINTER( ctypes.c_float ) )
	C23_ptr = C23.ctypes.data_as( ctypes.POINTER( ctypes.c_float ) )
	C33_ptr = C33.ctypes.data_as( ctypes.POINTER( ctypes.c_float ) )

	return C11_ptr,C12_ptr,C13_ptr,C22_ptr,C23_ptr,C33_ptr

def mngIdxs(idx,shp3):
	# Function to manage indexing of voxels to be processed in the 3D image (and associated Tensor).
	# This function is meant to harmonise indices values compatibly with a C-compiled shared library.
	#
	# Call: (idx,idx_ptr,numel,numel_ptr) = mngIdxs( idx,shp3 )
	#
	# *Inputs*
	# idx: scalar value identically equal to -1 (negative one)     *OR*
	#	   1D array containing unique C-like indices of selected voxel locations within the 3D image shape range.
	# shp3: shape of the 3D image
	#
	# *Outputs*
	# idx: updated 1D array containing unique C-like indices of 
	#	   *ALL* or *SELECTED* voxel locations within the 3D image shape range.
	# idx_ptr: C-like pointer to the 1D updated array of indices (idx)
	# numel: scalar number of elements in the updated 1D array of indices (idx)
	# numel_ptr: C-like pointer to the scalar number of elements (numel)

	# Import Libraries
	import numpy as np
	import ctypes
	
	# Checking Inputs
	if np.prod(np.shape(shp3)) != 3:
		raise Exception('Input: "shp3" should be the shape of a 3D image, e.g. (idx0, idx1, idx2). The value of "shp3" was: {}'.format(shp3))
	if ~ np.all( np.isfinite( shp3 ) , axis=-1 ):
		raise Exception('Input: "shp3" should contain all FINITE elements. The value of "shp3" was: {}'.format(shp3))
	if ~ np.all(shp3 == np.array(shp3,np.int32)):
		shp3 = np.array(shp3,np.int32)
		print(" -- Warning: Input argument (shp3) has non-integer elements - Default: round to closest integers.")
	if ~ np.all( np.greater(shp3 , 0) ):
		raise Exception('Input: "shp3" has non-positive number of elements. The value of "shp3" is: {}'.format(shp3))

	# Remove possible repetitions
	idx = np.unique(np.array(idx))

	if (np.prod(np.shape(idx)) == 1) & np.all(idx == -1):
		# CASE: process *ALL* voxels in the 3D image
		
		numel = np.array(np.prod(shp3), dtype=np.int32 ) #IMPORTANT: data type must be int32
		numel_ptr = numel.ctypes.data_as( ctypes.POINTER( ctypes.c_int ) )
		
		idx = np.array(np.arange(0,numel), dtype=np.int32 ) #IMPORTANT: data type must be int32
		idx_ptr = idx.ctypes.data_as( ctypes.POINTER( ctypes.c_int ) )
	
	else: 												  
		# CASE: process *SELECTED* voxels in the 3D image

		# Checking that the selected voxel indices are within the 3D image shape range
		if not ( np.all( np.greater_equal(idx , 0) ) & np.all( np.less(idx , np.prod(shp3)) ) ):
			raise Exception('The parsed indices (idx) are *OUT* of the 3D image. The parsed idx is: {}'.format( idx ) )

		idx = np.array(idx, dtype=np.int32 ) #IMPORTANT: data type must be int32
		idx_ptr = idx.ctypes.data_as( ctypes.POINTER( ctypes.c_int ) )
		
		numel = np.array(np.prod( idx.shape ), dtype=np.int32 ) #IMPORTANT: data type must be int32
		numel_ptr = numel.ctypes.data_as( ctypes.POINTER( ctypes.c_int ) )

	return idx,idx_ptr,numel,numel_ptr

def mpT3_eig(H11,H12,H13,H22,H23,H33,idx):
	# Function to determine the Eigen-Decomposition in 3D
	#
	# Call: (El1,El2,El3,Ev1,Ev2,Ev3,msk3Valid) = mpT3_eig( H11,H12,H13,H22,H23,H33,idx )
	#
	# *Inputs*
	# H11,H12,H13,H22,H23,H33 or [T11,T12,T13,T22,T23,T33]:
	# 	3D arrays with 6 independent (H)essian values or (T)ensorial components.
	# 	All H11,H12,...,H33 must have the same shape: shp3 = H11.shape, 
	# 	corresponding to a 3D image for each voxel (or location) at [i,j,k].
	# idx: 1-D array of C-like linear indices obtained from the function 'msk3_to_idx'
	#	   if idx is a scalar equal to -1, *ALL* indices are considered,
	#	   i.e. idx = np.array(np.arange(0,numel)), with numel = np.prod(shp3)
	#	   NB: voxel locations not selected by the array of indices (idx) will be considered
	#		   (and returned) as *isotropic* tensors by default.
	#
	# * Outputs *
	# El1,El2,El3: 3D arrays with values of shape shp3.
	# 	They represent the 3 Eigen-values associated to a
	# 	3D image at voxel (or location) [i,j,k].
	# Ev1,Ev2,Ev3: 4D arrays with values of shape (shp3,3).
	# 	They represent the 3 Eigen-vectors associated to the 
	# 	eigen-values of a 3D image at voxel (or location) [i,j,k].
	#
	# *** NB: Eigen-values and associated Eigen-vectors are sorted to satisfy
	#         the condition: abs(El1) <= abs(El2) <= abs(El3) 
	#
	# *** This is a python wrapper for the shared library: libmpT3libs.so ***
	# *** The C-compiled shared library automatically enables and configures
	# *** multi-core processing for maximal performance using OpenMP libraries.

	# Import Libraries
	import numpy as np
	import ctypes
	
	# Loading Shared Library
	mpT3libs = ctypes.cdll.LoadLibrary("./mpT3libs/libmpT3libs.so") #check path here!
	
	# Retrieving 3D image shape
	shp3 = np.shape(H11)
	# Managing voxel indices and pointers
	(idx,idx_ptr,
	 numel,numel_ptr) = mngIdxs(idx,shp3)

	# Regularising Inputs (variables and pointers)
	(H11,H12,H13,H22,H23,H33) = rgl6Cmps(H11,H12,H13,H22,H23,H33)
	# Determining Input Pointers
	(H11_ptr,H12_ptr,H13_ptr,
	 H22_ptr,H23_ptr,H33_ptr) = ptr6Cmps(H11,H12,H13,H22,H23,H33)

	# Initialising Outputs
	(El1,El2,El3,
	 Ev11,Ev12,Ev13,
	 Ev21,Ev22,Ev23,
	 Ev31,Ev32,Ev33) = iniElvs(shp3)
	# Determining Output Pointers 
	(El1_ptr,El2_ptr,El3_ptr,
	 Ev11_ptr,Ev12_ptr,Ev13_ptr,
	 Ev21_ptr,Ev22_ptr,Ev23_ptr,
	 Ev31_ptr,Ev32_ptr,Ev33_ptr) = ptrElvs(El1,El2,El3,
										   Ev11,Ev12,Ev13,
										   Ev21,Ev22,Ev23,
										   Ev31,Ev32,Ev33)

	# Initialising the Output Validity Mask
	(msk3Valid,msk3Valid_ptr) = inimsk3Valid(shp3)

	# Calling the C-compiled function from the shared library
	mpT3libs.mpT3_eig(H11_ptr,H12_ptr,H13_ptr,
					  H22_ptr,H23_ptr,H33_ptr,
					  El1_ptr,El2_ptr,El3_ptr,
					  Ev11_ptr,Ev12_ptr,Ev13_ptr,
					  Ev21_ptr,Ev22_ptr,Ev23_ptr,
					  Ev31_ptr,Ev32_ptr,Ev33_ptr,
					  msk3Valid_ptr, idx_ptr, numel_ptr)

	# Concatenating the Eigen-vectors components into 4D variables of shape (shp3,3)
	Ev1,Ev2,Ev3 = catEvs(Ev11,Ev12,Ev13,Ev21,Ev22,Ev23,Ev31,Ev32,Ev33)

	return El1,El2,El3,Ev1,Ev2,Ev3,msk3Valid

def mpT3_to_T3LIC(El1,El2,El3,Ev1,Ev2,Ev3,idx):
	# Function to convert the 3D Tensor field from its Eigen-Decomposition form
	# 	to its 6 (L)inearly (I)ndependent (C)omponent form.
	#
	# Call: (T11,T12,T13,T22,T23,T33) = mpT3_to_T3LIC( El1,El2,El3,Ev1,Ev2,Ev3,idx )
	#
	# *Inputs*
	# El1,El2,El3: 3D arrays with values of shape shp3.
	# 	They represent the 3 Eigen-values associated to a
	# 	Tensor field of a 3D image at voxel (or location) [i,j,k].
	# Ev1,Ev2,Ev3: 4D arrays with values of shape (shp3,3).
	# 	They represent the 3 Eigen-vectors associated to the respective
	# 	Eigen-values of a Tensor field of 3D image at voxel (or location) [i,j,k].
	# idx: 1-D array of C-like linear indices obtained from the function 'msk3_to_idx'
	#	   if idx is a scalar equal to -1, *ALL* indices are considered,
	#	   i.e. idx = np.array(np.arange(0,numel)), with numel = np.prod(shp3).
	#	   NB: voxel locations not selected by the array of indices (idx) will be considered
	#		   (and returned) as *isotropic* tensors by default.
	#
	# * Outputs *
	# T11,T12,T13,T22,T23,T33:
	# 	3D arrays with 6 independent (T)ensorial components.
	# 	All T11,T12,...,T33 must have the same shape (shp3), 
	# 	corresponding to a 3D image for each voxel (or location) at [i,j,k].
	# 
	# *** This is a python wrapper for the shared library: libmpT3libs.so ***
	# *** The C-compiled shared library automatically enables and configures
	# *** multi-core processing for maximal performance using OpenMP libraries.
	
	# Import Libraries
	import numpy as np
	import ctypes

	# Loading Shared Library
	mpT3libs = ctypes.cdll.LoadLibrary("./mpT3libs/libmpT3libs.so")

	# Retrieving 3D image shape
	shp3 = np.shape(El1)
	# Managing voxel indices and pointers
	(idx,idx_ptr,
	 numel,numel_ptr) = mngIdxs(idx,shp3)

	# Regularising inputs
	El1,El2,El3,Ev1,Ev2,Ev3 = rglElvs(El1,El2,El3,Ev1,Ev2,Ev3)
	# Separating Eigen-vectors components (for processing)
	(Ev11,Ev12,Ev13,
	 Ev21,Ev22,Ev23,
	 Ev31,Ev32,Ev33) = sepEvs(Ev1,Ev2,Ev3)
	# Determining Input Pointers
	(El1_ptr,El2_ptr,El3_ptr,
	 Ev11_ptr,Ev12_ptr,Ev13_ptr,
	 Ev21_ptr,Ev22_ptr,Ev23_ptr,
	 Ev31_ptr,Ev32_ptr,Ev33_ptr) = ptrElvs(El1,El2,El3,
										   Ev11,Ev12,Ev13,
										   Ev21,Ev22,Ev23,
										   Ev31,Ev32,Ev33)
	
	# Initialising outputs
	T11,T12,T13,T22,T23,T33 = ini6Cmps(shp3)
	# Determining Output Pointers 
	(T11_ptr,T12_ptr,T13_ptr,
	 T22_ptr,T23_ptr,T33_ptr) = ptr6Cmps(T11,T12,T13,T22,T23,T33)

	# Calling the C-compiled function from the shared library
	mpT3libs.mpT3_to_T3LIC(El1_ptr,El2_ptr,El3_ptr,
						   Ev11_ptr,Ev12_ptr,Ev13_ptr,
						   Ev21_ptr,Ev22_ptr,Ev23_ptr,
						   Ev31_ptr,Ev32_ptr,Ev33_ptr,
						   T11_ptr,T12_ptr,T13_ptr,
						   T22_ptr,T23_ptr,T33_ptr,
						   idx_ptr, numel_ptr)

	return T11,T12,T13,T22,T23,T33

def mpT3_to_T3LE(El1,El2,El3,Ev1,Ev2,Ev3,idx):
	# Function to convert the 3D Tensor field from its Eigen-Decomposition form
	# 	to its 6 linearly independent components in the (L)OG-(E)UCLIDEAN space.
	#
	# Call: (T11LE,T12LE,T13LE,T22LE,T23LE,T33LE) = mpT3_to_T3LE( El1,El2,El3,Ev1,Ev2,Ev3,idx )
	#
	# *Inputs*
	# El1,El2,El3: 3D arrays with values of shape shp3.
	# 	They represent the 3 Eigen-values associated to a
	# 	Tensor field of a 3D image at voxel (or location) [i,j,k].
	# Ev1,Ev2,Ev3: 4D arrays with values of shape (shp3,3).
	# 	They represent the 3 Eigen-vectors associated to the respective
	# 	Eigen-values of a Tensor field of 3D image at voxel (or location) [i,j,k].
	# idx: 1-D array of C-like linear indices obtained from the function 'msk3_to_idx'
	#	   if idx is a scalar equal to -1, *ALL* indices are considered,
	#	   i.e. idx = np.array(np.arange(0,numel)), with numel = np.prod(shp3).
	#
	# * Outputs *
	# T11LE,T12LE,T13LE,T22LE,T23LE,T33LE:
	# 	3D arrays with 6 independent (T)ensorial components in the LOG-EUCLIDEAN space.
	# 	All T11LE,T12LE,...,T33LE must have the same shape (shp3), 
	# 	corresponding to a 3D image for each voxel (or location) at [i,j,k].
	# 
	# *** This is a python wrapper for the shared library: libmpT3libs.so ***
	# *** The C-compiled shared library automatically enables and configures
	# *** multi-core processing for maximal performance using OpenMP libraries.

	# Import Libraries
	import numpy as np
	import ctypes

	# Loading Shared Library
	mpT3libs = ctypes.cdll.LoadLibrary("./mpT3libs/libmpT3libs.so")

	# Retrieving 3D image shape
	shp3 = np.shape(El1)
	# Managing voxel indices and pointers
	(idx,idx_ptr,
	 numel,numel_ptr) = mngIdxs(idx,shp3)

	# Regularising inputs
	El1,El2,El3,Ev1,Ev2,Ev3 = rglElvs(El1,El2,El3,Ev1,Ev2,Ev3)
	# Separating Eigen-vectors components (for processing)
	(Ev11,Ev12,Ev13,
	 Ev21,Ev22,Ev23,
	 Ev31,Ev32,Ev33) = sepEvs(Ev1,Ev2,Ev3)
	# Determining Input Pointers
	(El1_ptr,El2_ptr,El3_ptr,
	 Ev11_ptr,Ev12_ptr,Ev13_ptr,
	 Ev21_ptr,Ev22_ptr,Ev23_ptr,
	 Ev31_ptr,Ev32_ptr,Ev33_ptr) = ptrElvs(El1,El2,El3,
										   Ev11,Ev12,Ev13,
										   Ev21,Ev22,Ev23,
										   Ev31,Ev32,Ev33)

	# Initialising outputs
	T11LE,T12LE,T13LE,T22LE,T23LE,T33LE = ini6CmpsLE(shp3)
	# Determining Output Pointers 
	(T11LE_ptr,T12LE_ptr,T13LE_ptr,
	 T22LE_ptr,T23LE_ptr,T33LE_ptr) = ptr6Cmps(T11LE,T12LE,T13LE,T22LE,T23LE,T33LE)

	# Calling the C-compiled function from the shared library
	mpT3libs.mpT3_to_T3LE(El1_ptr,El2_ptr,El3_ptr,
						  Ev11_ptr,Ev12_ptr,Ev13_ptr,
						  Ev21_ptr,Ev22_ptr,Ev23_ptr,
						  Ev31_ptr,Ev32_ptr,Ev33_ptr,
						  T11LE_ptr,T12LE_ptr,T13LE_ptr,
						  T22LE_ptr,T23LE_ptr,T33LE_ptr,
						  idx_ptr, numel_ptr)

	return T11LE,T12LE,T13LE,T22LE,T23LE,T33LE

def mpT3LIC_to_T3(T11,T12,T13,T22,T23,T33,idx):
	# Function to convert the 3D Tensor field from its 6 (L)inearly (I)ndependent (C)omponent
	# form to its Eigen-Decomposition form, both in the EUCLIDEAN space.
	#
	# Call: (El1,El2,El3,Ev1,Ev2,Ev3,msk3Valid) = mpT3LIC_to_T3( T11,T12,T13,T22,T23,T33,idx )
	#
	# *Inputs*
	# T11,T12,T13,T22,T23,T33:
	# 	3D arrays with 6 linearly independent component of a (T)ensor field.
	# 	All T11,T12,...,T33 must have the same shape (shp3), 
	# 	corresponding to a 3D image for each voxel (or location) at [i,j,k].
	# idx: 1-D array of C-like linear indices obtained from the function 'msk3_to_idx'
	#	   if idx is a scalar equal to -1, *ALL* indices are considered,
	#	   i.e. idx = np.array(np.arange(0,numel)), with numel = np.prod(shp3).
	#	   NB: voxel locations not selected by the array of indices (idx) will be considered
	#		   (and returned) as *isotropic* tensors by default.
	#
	# * Outputs *
	# El1,El2,El3: 3D arrays with values of shape shp3.
	# 	They represent the 3 Eigen-values associated to a
	# 	Tensor field of a 3D image at voxel (or location) [i,j,k].
	# Ev1,Ev2,Ev3: 4D arrays with values of shape (shp3,3).
	# 	They represent the 3 Eigen-vectors associated to the respective
	# 	Eigen-values of a Tensor field of 3D image at voxel (or location) [i,j,k]. 
	# 
	# *** This is a python wrapper for the shared library: libmpT3libs.so ***
	# *** The C-compiled shared library automatically enables and configures
	# *** multi-core processing for maximal performance using OpenMP libraries.

	# Import Libraries
	import numpy as np
	import ctypes

	# Loading Shared Library
	mpT3libs = ctypes.cdll.LoadLibrary("./mpT3libs/libmpT3libs.so")

	# Retrieving 3D image shape
	shp3 = np.shape(T11)
	# Managing voxel indices and pointers
	(idx,idx_ptr,
	 numel,numel_ptr) = mngIdxs(idx,shp3)

	# Regularising Inputs
	T11,T12,T13,T22,T23,T33 = rgl6Cmps(T11,T12,T13,T22,T23,T33)
	# Determining Input Pointers
	(T11_ptr,T12_ptr,T13_ptr,
	 T22_ptr,T23_ptr,T33_ptr) = ptr6Cmps(T11,T12,T13,T22,T23,T33)

	# Initialising Outputs
	(El1,El2,El3,
	 Ev11,Ev12,Ev13,
	 Ev21,Ev22,Ev23,
	 Ev31,Ev32,Ev33) = iniElvs(shp3)
	# Determining Output Pointers 
	(El1_ptr,El2_ptr,El3_ptr,
	 Ev11_ptr,Ev12_ptr,Ev13_ptr,
	 Ev21_ptr,Ev22_ptr,Ev23_ptr,
	 Ev31_ptr,Ev32_ptr,Ev33_ptr) = ptrElvs(El1,El2,El3,
										   Ev11,Ev12,Ev13,
										   Ev21,Ev22,Ev23,
										   Ev31,Ev32,Ev33)

	# Initialising the Output Validity Mask
	(msk3Valid,msk3Valid_ptr) = inimsk3Valid(shp3)

	# Calling the C-compiled function from the shared library
	mpT3libs.mpT3LIC_to_T3(T11_ptr,T12_ptr,T13_ptr,
						   T22_ptr,T23_ptr,T33_ptr,
						   El1_ptr,El2_ptr,El3_ptr,
						   Ev11_ptr,Ev12_ptr,Ev13_ptr,
						   Ev21_ptr,Ev22_ptr,Ev23_ptr,
						   Ev31_ptr,Ev32_ptr,Ev33_ptr,
						   msk3Valid_ptr, idx_ptr, numel_ptr)

	# Concatenating the Eigen-vectors components into 4D variables of shape (shp3,3)
	Ev1,Ev2,Ev3 = catEvs(Ev11,Ev12,Ev13,Ev21,Ev22,Ev23,Ev31,Ev32,Ev33)

	return El1,El2,El3,Ev1,Ev2,Ev3,msk3Valid

def mpT3LE_to_T3(T11LE,T12LE,T13LE,T22LE,T23LE,T33LE,idx):
	# Function to convert the 3D Tensor field from its 6 linearly independent components
	# in the (L)OG-(E)UCLIDEAN space to its Eigen-Decomposition form in the EUCLIDEAN space.
	#
	# Call: (El1,El2,El3,Ev1,Ev2,Ev3) = mpT3LE_to_T3( T11LE,T12LE,T13LE,T22LE,T23LE,T33LE,idx )
	#
	# *Inputs*
	# T11LE,T12LE,T13LE,T22LE,T23LE,T33LE:
	# 	3D arrays with 6 independent (T)ensorial components in the LOG-EUCLIDEAN space.
	# 	All T11LE,T12LE,...,T33LE must have the same shape (shp3), 
	# 	corresponding to a 3D image for each voxel (or location) at [i,j,k].
	# idx: 1-D array of C-like linear indices obtained from the function 'msk3_to_idx'
	#	   if idx is a scalar equal to -1, *ALL* indices are considered,
	#	   i.e. idx = np.array(np.arange(0,numel)), with numel = np.prod(shp3).
	#	   NB: voxel locations not selected by the array of indices (idx) will be considered
	#		   (and returned) as *isotropic* tensors by default.
	#
	# * Outputs *
	# El1,El2,El3: 3D arrays with values of shape shp3.
	# 	They represent the 3 Eigen-values associated to a
	# 	Tensor field of a 3D image at voxel (or location) [i,j,k].
	# Ev1,Ev2,Ev3: 4D arrays with values of shape (shp3,3).
	# 	They represent the 3 Eigen-vectors associated to the respective
	# 	Eigen-values of a Tensor field of 3D image at voxel (or location) [i,j,k]. 
	# 
	# *** This is a python wrapper for the shared library: libmpT3libs.so ***
	# *** The C-compiled shared library automatically enables and configures
	# *** multi-core processing for maximal performance using OpenMP libraries.

	# Import Libraries
	import numpy as np
	import ctypes

	# Loading Shared Library
	mpT3libs = ctypes.cdll.LoadLibrary("./mpT3libs/libmpT3libs.so")

	# Retrieving 3D image shape
	shp3 = np.shape(T11LE)
	# Managing voxel indices and pointers
	(idx,idx_ptr,
	 numel,numel_ptr) = mngIdxs(idx,shp3)

	# Regularising Inputs
	T11LE,T12LE,T13LE,T22LE,T23LE,T33LE = rgl6Cmps(T11LE,T12LE,T13LE,T22LE,T23LE,T33LE)
	# Determining Input Pointers
	(T11LE_ptr,T12LE_ptr,T13LE_ptr,
	 T22LE_ptr,T23LE_ptr,T33LE_ptr) = ptr6Cmps(T11LE,T12LE,T13LE,T22LE,T23LE,T33LE)

	# Initialising Outputs
	(El1,El2,El3,
	 Ev11,Ev12,Ev13,
	 Ev21,Ev22,Ev23,
	 Ev31,Ev32,Ev33) = iniElvs(shp3)
	# Determining Output Pointers 
	(El1_ptr,El2_ptr,El3_ptr,
	 Ev11_ptr,Ev12_ptr,Ev13_ptr,
	 Ev21_ptr,Ev22_ptr,Ev23_ptr,
	 Ev31_ptr,Ev32_ptr,Ev33_ptr) = ptrElvs(El1,El2,El3,
										   Ev11,Ev12,Ev13,
										   Ev21,Ev22,Ev23,
										   Ev31,Ev32,Ev33)

	# Initialising the Output Validity Mask
	(msk3Valid,msk3Valid_ptr) = inimsk3Valid(shp3)

	# Calling the C-compiled function from the shared library
	mpT3libs.mpT3LE_to_T3(T11LE_ptr,T12LE_ptr,T13LE_ptr,
						  T22LE_ptr,T23LE_ptr,T33LE_ptr,
						  El1_ptr,El2_ptr,El3_ptr,
						  Ev11_ptr,Ev12_ptr,Ev13_ptr,
						  Ev21_ptr,Ev22_ptr,Ev23_ptr,
						  Ev31_ptr,Ev32_ptr,Ev33_ptr,
						  msk3Valid_ptr, idx_ptr, numel_ptr)

	# Concatenating the Eigen-vectors components into 4D variables of shape (shp3,3)
	Ev1,Ev2,Ev3 = catEvs(Ev11,Ev12,Ev13,Ev21,Ev22,Ev23,Ev31,Ev32,Ev33)

	return El1,El2,El3,Ev1,Ev2,Ev3,msk3Valid
