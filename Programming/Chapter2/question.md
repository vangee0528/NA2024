## 2.11 Problems

### 2.11.2 Programming assignments

A.Implement the Newton formula in a subroutine that produces the value of the interpolation polynomial $p_{n}(f;x_{0},x_{1},\ldots,x_{n};x)$ at any real x ,where $n\in\mathbb{N}^+$ ， $x_i'$ are distinct, and $f$ is a function assumed to be available in the form of a subroutine

B.Run your routine on the function

$$f(x)=\frac{1}{1+x^2}$$

for $x\in[-5,5]$ using $x_{i}=-5+10\frac{i}{n}$ ， $i=0,1,\ldots,n$ ,and $n$ = 2,4,6,8 .Plot the polynomials against the exact function to reproduce the plot in the notes that illus trate the Runge phenomenon

------------------------------------------------------------------

$$f(x)=\frac{1}{1+25x^2}$$

C. Reuse your subroutine of Newton interpolation to perform Chebyshev interpolation for the function

that extensively damage these trees in certain years.The following table lists the average weight of two samples of larvae at times in the first 28 days after birth. The first sample was reared on young oak leaves, whereas the sec ond sample was reared on mature leaves from the same tree.

for $x\in[-1,1]$ on the zeros of Chebyshev polynomials $T_{\mathrm{r}}$ with $n$ = 5, 10, 15, 20 .Clearly the Runge function $f(x)$ is a scaled version of the function in B. Plot the interpolating polynomials against the exact function to observe that the Chebyshev interpolation is free of the wide oscillations in the previous assignment

<table>
	<tbody>
		<tr>
			<th>$\overline{\mathrm{Day}}$</th>
			<th>0</th>
			<th>6</th>
			<th>$\overline{10}$</th>
			<th>$\overline{13}$</th>
			<th>$\overline{17}$</th>
			<th>$\overline{20}$</th>
			<th>28</th>
		</tr>
		<tr>
			<td>$\overline{\mathrm{Spl}}$</td>
			<td>6.67</td>
			<td>$\overline{17.3}$</td>
			<td>$\overline{42.7}$</td>
			<td>$\overline{37.3}$</td>
			<td>30.1</td>
			<td>29.3</td>
			<td>$\overline{28.7}$</td>
		</tr>
		<tr>
			<td>$\overline{\mathrm{Sp2}}$</td>
			<td>6.67</td>
			<td>$\overline{16.1}$</td>
			<td>18.9</td>
			<td>$\overline{15.0}$</td>
			<td>$\overline{10.6}$</td>
			<td>9.44</td>
			<td>8.89</td>
		</tr>
	</tbody>
</table>

D.A car traveling along a straight road is clocked at a num ber of points. The data from the observations are given in the following table, where the time is in seconds, the displacement is in feet, and the velocity is in feet per second.

(a) Use Newton's formula to approximate the average weight curve for each sample (b) Predict whether the two samples of larvae will die after another 15 days.

F.The roots of the following equation constitute a closed planar curve in the shape of a heart

<table>
	<tbody>
		<tr>
			<th>$\overline{\mathrm{Time}}$</th>
			<th>0</th>
			<th>3</th>
			<th>$\overline{5}$</th>
			<th>$\overline{8}$</th>
			<th>$\overline{13}$</th>
		</tr>
		<tr>
			<td>$\operatorname{displacement}$</td>
			<td>0</td>
			<td>$\overline{225}$</td>
			<td>383</td>
			<td>623</td>
			<td>$\overline{993}$</td>
		</tr>
		<tr>
			<td>velocity</td>
			<td>75</td>
			<td>777</td>
			<td>80</td>
			<td>74</td>
			<td>72</td>
		</tr>
	</tbody>
</table>

$$x^2+\left(\frac{3}{2}y-\sqrt{|x|}\right)^2=3.$$

(a) Use a Hermite polynomial to predict the position of the car and its speed for $t=10$ seconds

(b)Use the derivative of the Hermite polynomial to determine whether the car ever exceeds the speed limit of 55 miles per hour, i.e., 81 feet per second

Write a program to approximate the heart by cubic Bezier curves and plot your approximant, i.e. the piecewise cubic polynomials. Choose $m=10,40,160$ in Algo rithm 2.74 and produce three plots of the heart function Your knots should include the characteristic points

E.It is suspected that the high amounts of tannin in mature oak leaves inhibit the growth of the winter moth larvae