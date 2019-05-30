! STRUCT_HMS.INC

! Define some record structures

! ... where do we loose events?
	structure /STOP_EVENTS_H/
		integer*4	trials,slit_hor,slit_vert,slit_oct,
     >				Q1_in,Q1_mid,Q1_out,
     >				Q2_in,Q2_mid,Q2_out,
     >				Q3_in,Q3_mid,Q3_out,
     >				D1_in,D1_out,
     >				hut,incut,detec,successes
	end structure

! ... a general common block
	common /STOP_EVENTS/	STOP

	record /STOP_EVENTS_H/	STOP
