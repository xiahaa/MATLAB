function q = rot2quat(R)

	w = 0.5*sqrt( 1 + trace(R(1:3,1:3)) );
	w4 = w*4;
    q = [w;
		( R(3,2) - R(2,3) ) / w4;
		( R(1,3) - R(3,1) ) / w4;
        ( R(2,1) - R(1,2) ) / w4;
        ];
   
end