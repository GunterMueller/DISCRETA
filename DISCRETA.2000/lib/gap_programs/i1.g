# 
# i1.g
#
# Evi Haberberger
# December 1999
#


nu_prime := function(x, p)
	local i, q, r;
	
	i := 0;
	while x > 1 do
		r := RemInt(x, p);
		if r = 0 then
			i := i + 1; 
		else
			return i;
		fi;
		x := x/p;
	od;
	return i;
end;



find_indices := function(V, x)
# returns a vector of indices, where V has value x
	
	local l, i, res;
	
	res := [];
	l := Length(V);
	for i in [1..l] do
		if V[i] = x then
			Add(res, i);
		fi;
	od;
	return res;
end;

orbit_reps := function(G, W)
# returns orbit representatives w.r.t the action of G on W
	
	local W2, c, i, j, l, r, rep, p, perms, o, orbits, res;
	
	res := [];
	l := Length(W);
	W2 := [];
	perms := [];
	orbits := [];
	for i in [1..l] do
		Print("orbit_reps: step ", i, "\n");
		c := 0;
		for j in [1..i] do
		if c = 0 then
			# Print("conjugate? ", IsConjugate(G, W[i], W[j]), "\n");
			if i <> j and IsConjugate(G, W[i], W[j]) = true then
				c := c + 1;
			fi;
		fi;
		od;
		if c = 0 then
			p := [()];
			o := [W[i]];
			Add(W2, W[i]);
			# Print("representative: ", W[i], "\n");
			for j in [i+1..l] do
				if IsConjugate(G, W[i], W[j]) = true then
					Add(o, W[j]);
					Add(p, RepresentativeOperation(G, W[i], W[j]));
				fi;
			od;
			Add(orbits, o);
			# Print("orbits[", i, "]:", o);
			Add(perms, p);
			# Print("permutations[", i, "]:", p, "\n");
		fi;
	od;
	Add(res, W2);
	Add(res, orbits);
	Add(res, perms);
	return res;
end;

is_constant := function( V, S )
# checks if the vector V is constant on the set of indices S,
# returns -1 if not, or the value c which V has on all these positions.

	local i, idx, l, c;

	l := Length(S);
	if l = 0 then
  		Print("ERROR: is_constant, indexset is empty\n");
	fi;
	for i in [1..l] do
  		idx := S[i];
  		if i = 1 then
    			c := V[idx];
  		else
  	  		if V[idx] <> c then
      				return -1;
    			fi;
  		fi;
	od;
	return c;
end;


create_overgroups := function(r, Agen, G, p, n0, strtex, fv)
	local z, l, h, Agenh, lA, j, a, Ah, B, b, n, f, W, ind, grsize, res;
	
	Print("\nstart create_overgroups():\n\n");
	res := [];
	l := Length(r);
	W := [];
	ind := [];
	grsize := [];
	for z in [1..l] do
		h := r[z];
		# Print("h := ", h, "\n");
		Agenh := [];
		lA := Length(Agen);
		for j in [1..lA] do
			a := Agen[j]^h;
			Add(Agenh, a);
		od;
		Ah := Concatenation(Agen, Agenh);
		B := Subgroup(G, Ah);
		b := Size(B);
		# Print("group order ", z, ": ", b, "\n");
		n := nu_prime(b, p);
		f := Size(G)/b;
		if f > 2 then
			if IsBound(W[n]) = false then
				W[n] := [];
				grsize[n] := [];
				ind[n] := [];
			fi;
			Add(W[n], B);
			Add(grsize[n], b);
			Add(ind[n], z);
			Print(z, " ", n, "\n");
		fi;
	od;
	res[1] := Compacted(W);
	res[2] := Compacted(grsize);
	res[3] := Compacted(ind);
	Print("\nleave create_overgroups()\n\n");
	return res;
end;
	
action_by_conjugation := function(NGP, W, strtex, fv)
# action of N_G(P) on the B by conjugation (this is possible,
# when all of them have P as a common p-Sylowgroup!):

	local rep, ll, res, Brep, Borb, Bperm, x, lll, Bperm1, Borb1, 
	eq, y, z, i, li;
	
	rep := orbit_reps(NGP, W);
	# rep[1]: representatives, 
	# rep[2]: conjugates of rep[1],
	# rep[3]: conjugation permutation (first is identity)
	
	Print("\nstart action_by_conjugation(): orbit reps\n", W, "\n");
	
	ll := Length(rep[1]);
	Print("number of conjugacy classes under NGP:\n", ll, "\n\n");
	res := [];
	Brep := [];
	Borb := [];
	Bperm := [];
	for x in [1..ll] do
		lll := Length(rep[3][x]);
		Add(Brep, rep[1][x]);
		Bperm1 := [];
		Borb1 := [];
		for z in [1..lll] do
			eq := 0;
			for y in [1..z] do
			if y < z then
				if rep[2][x][y] = rep[2][x][z] then
					eq := eq + 1;
				fi;
			fi;
			od;
			# Print("IsZero?", eq, "\n");
			if eq = 0 then
				Add(Borb1, rep[2][x][z]);
				Add(Bperm1, rep[3][x][z]);
			fi;
		od;
		# Print("Borb = ", Borb, "\n");
		Add(Borb, Borb1);
		Add(Bperm, Bperm1);
	od;

	# Print("Brep = ", Brep, "\n");
	ll := Length(Brep);
	if fv > 0 then
		AppendTo(strtex, "\n\n\\medskip\nWe receive a partition of those groups into ");
		AppendTo(strtex, String(ll));
		AppendTo(strtex, " conjugacy classes \n(the first group is $A$ itself).\n");
		if fv > 1 and ll > 0 then
			AppendTo(strtex, "\n\\begin{itemize}\n");
			for i in [1..ll] do
				AppendTo(strtex, "\n\\item\n");
				AppendTo(strtex, "representative ");
				AppendTo(strtex, String(i));
				AppendTo(strtex, ": ");
#				AppendTo(strtex, Brep[i]);
				if IsSolvableGroup(Brep[i]) = 0 then
					AppendTo(strtex, "solvable ");
				# elif IsSymmetricGroup(Brep[i]) = 0 then
				#	AppendTo(strtex, "symmetric ");
				# elif IsAlternatingGroup(Brep[i]) = 0 then
				# 	AppendTo(strtex, "alternating ");
				elif IsPerfectGroup(Brep[i]) = 0 then
					AppendTo(strtex, "perfect ");
				fi;
				AppendTo(strtex, "\\\\ \ngroup order: ");
				AppendTo(strtex, String(Size(Brep[i])));	
				AppendTo(strtex, "\\\\ \nnumber of conjugates (including the representative): ");
				li := Length(Borb[i]);
				AppendTo(strtex, String(li));
			od; 
			AppendTo(strtex, "\n\\end{itemize}\n");
		fi;
		AppendTo(strtex, "\n\\medskip\n");
	fi;

	Add(res, Brep);
	Add(res, Borb);
	Add(res, Bperm);
	Print("leave action_by_conjugation()\n\n");
	return res;
end;
