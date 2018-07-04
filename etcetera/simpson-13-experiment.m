format long;
%{
--	First item:
	"Escreva um programa na linguagem C ou matlab que
	utiliza o metodo de Simpson 1/3 composta para cal-
	cular um valor aproximado da integral"
%}

%{
	Problem formula
%}
function res = f (x)
	res = exp(-x**2 / 2)/(2 * pi)**0.5;
endfunction

%{
	Simpson 1/3 formula:
	integrate{a <= x <= b}f(x)dx approx 
		(b-a)/6 * (f(a) + 4*f((a + b)/2) + f(b))
%}

function I = simpson (a = 0, b = 1, fun = @f)
	I = (b-a)/6 * (fun(a) + 4*fun((a+b)/2) + fun(b));
endfunction

%{
	F(x) = I(f)(x) - 0.45
%}
function res = F(a = 0, b = 1, fun = @f) 
	res = simpson(a, b, fun) - 0.45;
endfunction

%{
--	Second item:
	Use o programa acima e verifique que F(1)F(2)<0, em que 
	F(1) = I(f)(1)−0.45, idem para F(2).  Isso determina o 
	intervalo onde F(x) possui raiz.
%}

F(a = 0, b = 1)*F(a = 0, b = 2) < 0 
% True!

%{
-- Final item:
	Utilize o metodo de Newton e determine a raiz z in [1,2] 
	da Eq. 1 com precisao EPS = 10^{−10} (tome x_{0} = 0.5 
	e CALCULE ITERACOES x_{n+1} ate que |x_{n+1}−x_n|< EPS).
%}

%{
	Newton formula:
	x_{n+1} = x_{n} - F(x)/f(x)	,	f(x) = F'(x)

%}
function res = newton (x0 = .5, EPS = 10e-10, fun = @f)
	err = 2.0*EPS;
	x_cur = x0;

	while (err > EPS)
		x_new = x_cur - F(a = 0, b = x_cur, fun = fun)/f(x_cur);
		err = abs(x_new - x_cur)
		x_cur = x_new
	endwhile

	res = x_cur;
endfunction

res = newton()

% "Real proof"
simpson(a=0, b=res) % = .45
