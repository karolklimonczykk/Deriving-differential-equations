# Deriving differential equations
Deriving differential equations of motion of discrete systems (e.g. double pendulum). To derive the equations of motion, I'm using the Lagrangian procedure, which generates differential equations of motion based on Lagrangian equations of the second kind.<br><br>

<code>lagrange := proc(n, q, r, L)</code> <br>
&nbsp; &nbsp;<code>local i, uzm_q, uzm_r, rel_r_q, Lq, Lr, Lrt;</code><br>
&nbsp; &nbsp;<code>global row;</code><br>
 &nbsp; &nbsp;<code>uzm_q := seq(q[i] = q[i] (t), i = 1 .. n);</code><br>
  &nbsp; &nbsp;<code>uzm_r := seq(r[i] = r[i] (t), i = 1 .. n);</code><br>
  <i>differentiation by generalized coordinates and their variation</i><br>
 &nbsp; &nbsp;<code>for i to n do</code><br>
 &nbsp; &nbsp;&nbsp; &nbsp;<code>Lq[i] := subs([uzm_q, uzm_r], diff(L, q[i]));</code><br>
&nbsp; &nbsp;&nbsp; &nbsp;<code>Lr[i] := subs([uzm_q, uzm_r], diff(L, r[i]));</code><br>
&nbsp; &nbsp;<code>end do;</code> <br>
<i>differentiation over time of the derivatives of L with respect to velocities</i><br>
  &nbsp; &nbsp;<code>for i to n do</code><br>
   &nbsp; &nbsp;&nbsp; &nbsp;<code>Lrt[i] := diff(Lr[i], t);</code><br>
&nbsp; &nbsp;<code>end do;</code><br> 
<i>substitution of mutual relations between velocities and generalized displacements</i><br>
 &nbsp; &nbsp;<code>rel_r_q := seq(r[i] (t) = diff(q[i] (t), t), i = 1 .. n);</code><br>
 <i>generating the equation of motion</i><br>
  &nbsp; &nbsp;<code>for i to n do</code><br>
   &nbsp; &nbsp;&nbsp; &nbsp;<code>row[i] := subs(rel_r_q, Lrt[i] - Lq[i] = 0);</code><br>
  &nbsp; &nbsp;<code>end do;</code><br>
  &nbsp; &nbsp;<code>seq(row[i], i = 1 .. n);</code><br>
<code>end proc</code>
<br><br>
<p>The formal parameters of the procedure are:</p>
<ul>
<li>n – number of degrees of freedom,</li>
 <li>q – name of the indexed variable defining generalized coordinates,</li>
 <li>r – name of the indexed variable defining generalized velocities,</li>
 <li>L – Lagrange function.</li> <br>
 </ul>
 The global variables row[1], ..., row[n] were assigned the searched equations of motion.
