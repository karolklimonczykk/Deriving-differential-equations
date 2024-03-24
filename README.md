# Deriving differential equations
Deriving differential equations of motion of discrete systems (e.g. double pendulum). To derive the equations of motion, I'm using the Lagrangian procedure, which generates differential equations of motion based on Lagrangian equations of the second type.
<code>
lagrange := proc(n, q, r, L)</code> <br>
&nbsp; &nbsp;<code>local i, uzm_q, uzm_r, rel_r_q, Lq, Lr, Lrt;</code><br>
&nbsp; &nbsp;<code>global row;</code><br>
 &nbsp; &nbsp;<code>uzm_q := seq(q[i] = q[i](t), i = 1 .. n);</code><br>
  &nbsp; &nbsp;<code>uzm_r := seq(r[i] = r[i](t), i = 1 .. n);</code><br>
 &nbsp; &nbsp;<code>for i to n do</code><br>
 &nbsp; &nbsp;&nbsp; &nbsp;<code>Lq[i] := subs([uzm_q, uzm_r], diff(L, q[i]));</code><br>
&nbsp; &nbsp;&nbsp; &nbsp;<code>Lr[i] := subs([uzm_q, uzm_r], diff(L, r[i]));</code><br>
&nbsp; &nbsp;<code>end do;</code> <br>
  &nbsp; &nbsp;<code>for i to n do</code><br>
   &nbsp; &nbsp;&nbsp; &nbsp;<code>Lrt[i] := diff(Lr[i], t);</code><br>
&nbsp; &nbsp;<code>end do;</code><br> 
 &nbsp; &nbsp;<code>rel_r_q := seq(r[i](t) = diff(q[i](t), t), i = 1 .. n);</code><br>
  &nbsp; &nbsp;<code>for i to n do</code><br>
   &nbsp; &nbsp;&nbsp; &nbsp;<code>row[i] := subs(rel_r_q, Lrt[i] - Lq[i] = 0);</code><br>
  &nbsp; &nbsp;<code>end do;</code><br>
  &nbsp; &nbsp;<code>seq(row[i], i = 1 .. n);</code><br>
<code>end proc</code>
