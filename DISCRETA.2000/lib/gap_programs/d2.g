# d2.g
#
# Evi Haberberger
# August/September 1999
#


topological_sort := function(Bsol, Bcon, param)
	
	local L, l, i, z, B, b, ol, sl, info, grlist, vec, lvec, top,
		b0, vec0, lvec0, j, bj, vecj, lvecj, ll, k, test, res,
		len, len1, g, f;
	
	L := Bsol[1];
	l := Length(Bsol[1]);
	res := [];
	grlist := [];
	for i in [1..l] do
		info := [];
		z := L[i];		# KM-identifier
		Add(info, z);
		B := Bcon[1][z];	# conjugacy representative (order saved)
		b := Size(B);
		Add(info, b);
		ol := Length(Bcon[2][z]); # nb. of conjugates including rep.
		Add(info, ol);
		
		sl := Length(Bsol[3][i]); # nb. of solutions for each conj.
		Add(info, sl);
		Add(grlist, info);
	od;
	
	# now sort topologically:
	top := [];
	b0 := grlist[1][2];
	vec0 := FactorsInt(b0);
	lvec0 := Length(vec0);
	top[lvec0] := [];
	Add(top[lvec0], grlist[1]);
	for i in [2..l] do
		b := grlist[i][2];
		vec := FactorsInt(b);
		lvec := Length(vec);
		if IsBound(top[lvec]) = true then
			ll := Length(top[lvec]);	
			# top[lvec] is lattice layer
		else
			ll := 0;
			top[lvec] := [];
		fi;
		test := 0;
		for j in [1..ll] do
			bj := top[lvec][j][2];
			if bj > b then
				for k in [j..ll] do
					top[lvec][ll+j-k+1] := top[lvec][ll+j-k];
				od;
				top[lvec][j] := grlist[i];
				test := test+1;
			fi;
		od;
		if test=0 then	# not yet added (->to the end)! 
			if ll = 0 then
				top[lvec] := [];
			fi;
			Add(top[lvec], grlist[i]);
		fi;
	od;
	res := Compacted(top);
	Print("res=", res, "\n");
	
	# in the end add the symmetric group G itself:
	len := Length(res);
	Add(res, []);
	g := [];
	l := Length(Bcon[1]);
	Add(g, l+1);
	f := Factorial(param[1]); 
	Add(g, f);	# group order of S_v
	Add(g, 1);	# no other conjugates
	Add(g, 0);	# no solutions
	Add(res[len+1], g);
	
	# add the layer number for each group
	len := Length(res);
	for i in [1..len] do
		len1 := Length(res[i]);
		for j in [1..len1] do
			Add(res[i][j], i);
		od;
	od;
	
	return res;
end;

create_Ainf := function(top, Bcon, param, strmp)
	
	local Ainf, len, l, i, ll, j, z, B, e, ii, jj, lll, zz, ol, k, BB, x, 
	y;
	
	Ainf := [];
	l := Length(top);
	len := 0;
	# add the group infos to the mpbase-file:
	for i in [1..l] do
		Print(top[i]);
		ll := top[i][1][3];
		len := len + ll;
	od;
	AppendTo(strmp, String(l));
	AppendTo(strmp, " ");
	AppendTo(strmp, String(len));
	AppendTo(strmp, "\n\n");
	for i in [1..l] do
		ll := Length(top[i]);
		for j in [1..ll] do
			for k in [1..5] do;
				AppendTo(strmp, String(top[i][j][k]));
				AppendTo(strmp, " ");
			od;
			AppendTo(strmp, "\n");
		od;
	od;
	AppendTo(strmp, "\n");
	# initialize Ainf:	
	for i in [1..l] do
		Ainf[i] := [];
		for j in [1..l] do
			if i<>j then
				Ainf[i][j] := 0;
			else
				Ainf[i][j] := 1;
			fi;
		od;
	od;
	
	# now start the computation:
	x := 0;
	for i in [1..l-1] do
		Print("layer=", i, "\n");
		ll := Length(top[i]);
		# Print("length=", ll, "\n");
		for j in [1..ll] do
			# Print("group ", j, ":\n");
			z := top[i][j][1];
			B := Bcon[1][z];
			# test the elements of the upper conjugacy classes,
			# if they are overgroups:
			y := x;
			for ii in [i+1..l] do
				Print("tested layer:", ii, "\n");
				if ii < l then
					lll := Length(top[ii]);
					# Print("length=", lll, "\n");
					for jj in [1..lll] do
						# Print("group ", jj, ":\n");
						zz := top[ii][jj][1];
						e := 0;
						# ol = length of conj. class
						ol := Length(Bcon[2][zz]);
						for k in [1..ol] do
							BB := Bcon[2][zz][k];
							# Print("group ", k, "\n");
							if IsSubgroup(BB, B)=true then
								e := e+1;
							fi;
						od;
						# Print("row=", x+j, " column=", y+jj+1, ": entry=",e, "\n"); 
						Ainf[x+j][y+jj+1] := e;
					od;
				else
					lll := 1;
					BB := SymmetricGroup(param[1]);
					if IsSubgroup(BB, B)=true then
						Ainf[x+j][y+lll+1] := 1;
					fi; 
				fi;
				y := y+lll;
			od;
		od;
		x := x+ll;
	od;
	Print("Ainf=", Ainf, "\n");
	for i in [1..l] do
		for j in [1..l] do
			AppendTo(strmp, String(Ainf[i][j]));
			AppendTo(strmp, " ");
		od;
		AppendTo(strmp, "\n");
	od;
	AppendTo(strmp, "\n");
	return Ainf;
end;

plesken_submatrix_ij := function(top, Bcon, param, m, n, strmp)
	
	local pl, l, ol1, ol2, i, j, zi, zj, ki, kj, Bi, Bj;
	
	pl := [];
	pl[1] := m;	# submatrix of plesken: for (m, n) with m < n
	pl[2] := n;
	
	# we test the subgroup relations between 
	# conjugacy class m and n (if they are in different layers)
	
	AppendTo(strmp, String(m));
	AppendTo(strmp, " ");
	AppendTo(strmp, String(n));
	AppendTo(strmp, "\n");
	
	Add(pl, []);	
	if top[m][5] = top[n][5] then	# are in the same layer -> 0-matrix
		ol1 := top[m][3];
		ol2 := top[n][3];
		for i in [1..ol1] do
			Add(pl[3], []);
			for j in [1..ol2] do
				Add(pl[3][i], 0);
				AppendTo(strmp, 0);
				AppendTo(strmp, " ");
			od;
			AppendTo(strmp, "\n");
		od;
	else
		if n < Length(top) then
			zi := top[m][1];
			zj := top[n][1];
			ol1 := Length(Bcon[2][zi]);
			ol2 := Length(Bcon[2][zj]);
			for ki in [1..ol1] do
				Bi := Bcon[2][zi][ki];
				Add(pl[3], []);
				for kj in [1..ol2] do
					Bj := Bcon[2][zj][kj];
					if IsSubgroup(Bj, Bi)=true then
						Add(pl[3][ki], 1);
						AppendTo(strmp, 1);
						AppendTo(strmp, " ");
					else
						Add(pl[3][ki], 0);
						AppendTo(strmp, 0);
						AppendTo(strmp, " ");
					fi;
				od;
				AppendTo(strmp, "\n");
			od; 
		else
			zi := top[m][1];
			ol1 := Length(Bcon[2][zi]);
			for ki in [1..ol1] do
				Bi := Bcon[2][zi][ki];
				Add(pl[3], []);
				Bj := SymmetricGroup(param[1]);
				if IsSubgroup(Bj, Bi)=true then
					Add(pl[3][ki], 1);
					AppendTo(strmp, 1);
					AppendTo(strmp, " ");
				else
					Add(pl[3][ki], 0);
					AppendTo(strmp, 0);
					AppendTo(strmp, " ");
				fi;
				AppendTo(strmp, "\n");
			od;
		fi;
	fi;
	AppendTo(strmp, "\n");
	return pl;
end;

group_sol_lattice := function(Bsol, Bcon, param, strmp)
	
	local l, lt, i, top, Ainf, ll, top0, j, lll, len;
	
	Print("\nstart group_sol_lattice():\n\n");
	l := Length(Bsol[1]);	# nb of group reps with solutions
	top := topological_sort(Bsol, Bcon, param);
	Ainf := create_Ainf(top, Bcon, param, strmp);
	
	# first modify top: eliminate the partitioning into layers
	ll := Length(top);
	top0 := [];
	for i in [1..ll] do
		lll := Length(top[i]);
		for j in [1..lll] do
			Add(top0, top[i][j]);
		od;
	od;
	len := Length(top0);
	for i in [1..len] do
		for j in [i+1..len] do
			plesken_submatrix_ij(top0, Bcon, param, i, j, strmp);
		od;
	od;
	Print("\nleave group_sol_lattice()\n\n");
	return top;
end;
