# 
# i3.g
#
# Evi Haberberger
# December 1999
#



action_N_NGP_A_2 := function(N_NGP_A, canon, param, KM, strtex, fv) 
# Q gives the operation of N_NGP_A on the orbits and hence the resulting solutions:

	local lg, Q, q, qq, q1, km, km1, i, str, l, j; 

	Print("\n\nstart action_N_NGP_A_2():\n");
	lg := Length(GeneratorsOfGroup(N_NGP_A));
	if fv > 0 then
		AppendTo(strtex, "\\subsection{Orbits of the action of $N_{NGP}(A)$ on the solutions}\n\n");
		if fv > 1 then
			AppendTo(strtex, "The action of $N_{NGP}(A)$ on the solutions, which are ");
			AppendTo(strtex, "only invariant under $A$, \nis determined by the ");
			AppendTo(strtex, "following permutations in $Q$:\n\n\\smallskip\n");
			AppendTo(strtex, "%{\\tiny \n");
		fi;
	fi;
	Q := [];
	qq := PermList(canon[2][1]);
	for i in [1..lg] do
	  	km := [];
		Append(km, String(KM));
		q := action_on_blocks(km, param[3], GeneratorsOfGroup(N_NGP_A)[i]);
		q1 := qq^-1 * q * qq;
		Add(Q, q1);
		str := String(q1);
		if fv > 2 then
			AppendTo(strtex, "\\begin{align*}\n");
			AppendTo(strtex, str);
			AppendTo(strtex, "\n\\end{align*}\n");
		fi;
	od;
	if fv > 2 then
		AppendTo(strtex, "\n%}\n\n\\medskip\n");
	fi;
	Print("\nleave action_N_NGP_A_2()\n\n");
	return Q;
end;
