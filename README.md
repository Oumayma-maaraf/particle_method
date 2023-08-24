\textbf{Stochastic collocation}

Soit Y une variable aléatoire à valeurs réelles avec une fonction de répartition strictement croissante F_{Y}(y). Soit U\sim\mathcal{U}([0,1]) et u_{n}un échantillon de U. Classiquement pour générer un échantillon de Y on fait 

y_{n}=F_{Y}^{-1}(u_{n})

Or dans le cas où la fonction inverse de répartition n'a pas de forme analytique cette procédure devient couteuse.

On alors considère une autre variable X, pour laquelle F_{X}^{-1}(u_{n}) est moins couteuse que celle de Y. On sait que F_{Y}(Y)\overset{d}{=}F_{X}(X) donc y_{n}=F_{Y}^{-1}(F_{X}(\xi_{n})) où y_{n},\xi_{n} les échantillons de Y,X respectivement. Ici encore l'échantillonement de Y reste couteux. Il faut donc trouver une relation alternative pour ne pas faire l'invesion F_{Y}^{-1}pour tous l'échantillonf de X.

On cherche donc une fonction de manière à ce que g(.)=F_{Y}^{-1}(F_{X}(.)) donc à ce que Y\overset{d}{=}g(X), et telle que l'évaluation de cette fonction n'est pas couteuse.

Dans la méthode de collocation stochastique on approxime Y par une foncition g de X en terme d'expansion de Lagrange l_{i}(\xi_{n}):y_{n}\approx g_{N}(\xi_{n})=\sum_{i=1}^{N}y_{i}l_{i}(\xi_{n}),\quad l_{i}(\xi_{n})=\prod_{j=1,i\neq j}^{N}\frac{\xi_{n}-x_{j}}{x_{i}-x_{j}}

où \xi_{n}est un échantillon de X et x_{i},x_{j}sont des points de collocation, et y_{i}=F_{Y}^{-1}(F_{X}(x_{i})). \textbf{l}(x)=(l_{1}(x),\dots,l_{N}(x))^{T} est la base de Lagrange., telle que l_{i}(x_{j})=\delta_{ij}. Donc une fois les N points de collocation déterminés x_{i} et les N inversions F_{Y}^{-1} faites, on peut simuler nimporte quel nombre d'échantillons de Y et ceci par l'évaluation du polynome g_{N}(.). On parle ici de \textit{Stochastic Collocation Monte Carlo sampler}.

• \textbf{Points de collocation}

On sait que g_{N}(x)=\sum_{i=1}^{N}y_{i}l_{i}(x),\quad l_{i}(x)=\prod_{j=1,i\neq j}^{N}\frac{x-x_{j}}{x_{i}-x_{j}}

Poue éviter de faire O(N^{2}) opération pour chaque nouvelle valeurs de x, on considère un poids \lambda_{i} défini par 

\lambda_{i}=\frac{1}{\prod_{j=1,j\neq i}^{N}(x_{i}-x_{j})},\quad i\in\llbracket1,N\rrbracket

et l(x)=(x-x_{1})\dots(x-x_{N}), donc l_{i}(x) peut s'écrire comme

l_{i}(x)=l(x)\frac{\lambda_{i}}{x-x_{i}}

donc g_{N}(x)=\sum_{i=1}^{N}\frac{y_{i}\lambda_{i}}{(x-x_{i})}l(x)

On dit qu'une sequence de polynomes orthogonaux \{p_{i}\}_{i=0}^{N} avec degré deg(p_{i})=i est orthogonale en L^{2}par rapport à la densité de f_{X}(X) de X, si

\mathbb{E}(p_{i}(X)p_{j}(X))=\delta_{ij}\mathbb{E}(p_{i}^{2}(X)),\quad i,j=0,\dots,N

Et pour toutes densité f_{X}(.), il existe une sequence de polynomes orthogonaux p_{i}(x) unique avec comme degré deg(p_{i}(x))=i, cette sequence se construit par 

p_{i+1}(x)=(x-\alpha_{i})p_{i}(x)-\beta_{i}p_{i-1}(x),\quad i\in\llbracket0,N-1\rrbracket

où p_{-1}(x)=0 et p_{0}(x)=1 et pour tous i\in\llbracket0,N-1\rrbracket \alpha_{i}=\frac{\mathbb{E}(Xp_{i}^{2}(X))}{\mathbb{E}(p_{i}^{2}(X))} et \beta_{i}=\frac{\mathbb{E}(p_{i}^{2}(X))}{\mathbb{E}(p_{i-1}^{2}(X))} et \beta_{0}=0

On construit la matrice de Gram M=\{\mu_{ij}\}_{i,j=0}^{N} en considérant le monomial m_{i}(X)=X^{i} par \mu_{ij}=\mathbb{E}(m_{i}(X)m_{j}(X))=\mathbb{E}(X^{i+j}). La matrice M est définie positive elle s'écrit alors comme M=R^{T}R où R=\{r_{ij}\}_{i,j=0}^{N} est une matrice triangulaire inférieure. On a 

\alpha_{j}=\frac{r_{j,j+1}}{r_{j,j}}-\frac{r_{j-1,j}}{r_{j-1,j-1}},\quad\beta_{j}=\left(\frac{r_{j+1,j+1}}{r_{j,j}}\right)^{2},\quad j=1,\dots,N-1

où r_{0,0}=1 et r_{0,1}=0

A.2- Les zeros x_{i},i\in\llbracket1,N\rrbracket du polynôme orthogonal p_{N}(X) sont les valeurs propres de la matrice symétrique suivante

\hat{J}:=\left(\begin{array}{ccccc}
\alpha_{1} & \sqrt{\beta_{1}} & 0 & 0 & 0\\
\sqrt{\beta_{1}} & \alpha_{2} & \sqrt{\beta_{2}} & 0 & 0\\
 & \ddots & \ddots & \ddots\\
0 & 0 & \sqrt{\beta_{N-2}} & \alpha_{N-1} & \sqrt{\beta_{N-1}}\\
0 & 0 & 0 & \sqrt{\beta_{N-1}} & \alpha_{N}
\end{array}\right)

• \textbf{Analyse d'erreur}

On s'intéresse dans cette section à l'erreur générée par \textit{Stochastic Collocation Monte Carlo sampler}.

On se met dans un premier temps dans un cas où la méthode de collocation donne des résultats exacts. Soient Y\sim\mathcal{N}(\mu_{Y},\sigma_{Y}^{2}) et X\sim\mathcal{N}(\mu_{X},\sigma_{X}^{2}) deux variables aléatoires alors g_{N}(X)\overset{d}{=}Y pour N=2. En effet, soietn x_{1}et x_{2} deux points de collocation alorsg_{2}(X)=y_{1}\frac{X-x_{2}}{x_{1}-x_{2}}+y_{2}\frac{X-x_{1}}{x_{2}-x_{1}}

On a F_{\mathcal{N}(0,1)}\left(\frac{y_{i}-\mu_{Y}}{\sigma_{Y}}\right)=F_{\mathcal{N}(0,1)}\left(\frac{x_{i}-\mu_{X}}{\sigma_{X}}\right), alors y_{i}=\frac{x_{i}-\mu_{X}}{\sigma_{X}}\sigma_{Y}+\mu_{Y}. On a donc \mathbb{E}\left(g_{2}(X)\right)=\mu_{Y} et \mathbb{V}\left(g_{2}(X)\right)=\sigma_{Y}^{2} et comme g_{2}(X) suit une loi normale alors Y\overset{d}{=}g_{2}(X).

Dans un cas général, pour mesurer l'erreur on peut soit considérer la différence entre g(X) et g_{N}(X) soit l'erreur associée à l'approximation de la fonction de répartition.

La première erreur est liée à l'interpolation de Lagrange, en effet la relation entre Y et X est Y=g(X) qu'on approxime par un polynome de Lagrange Y\approx g_{N}(X) pour N points de collocation. Cette erreur est donc bien connue 

e_{X}(\xi_{n})=|g(\xi_{n})-g_{N}(\xi_{n})|=\bigg|\frac{1}{N!}\frac{d^{N}g(x)}{dx^{N}}\bigg|_{x=\hat{\xi}}\prod_{i=1}^{N}(\xi_{n}-x_{i})\bigg|

avec x_{i}est un point de collocation, \hat{\xi}\in[\min(x),\max(x)] et x=(x_{1},\dots,x_{N})^{T}, on peut borner cette erreur en prenant \hat{\xi} l'abscisse du maximum de \bigg|\frac{d^{N}g(x)}{dx^{N}}\bigg|. 

\textbf{Erreur de convergence en}L^{2}

On a Y=g(X)\approx Y_{N}\equiv g_{N}(X), où g(x)=F_{Y}^{-1}(F_{X}(x)) donc 

\mathbb{E}\left[(Y-Y_{N})^{2}\right]=\mathbb{E}\left[(g(X)-g_{N}(X))^{2}\right]=\int_{\mathbb{R}}(g(x)-g_{N}(x))^{2}f_{X}(x)dx

Les points de collocations x_{i} et les poinds w_{i}sont déterminés par le théorème A.2. Comme g(x_{i})=g_{N}(x_{i}),pour i\in\llbracket1,N\rrbracket l'erreur est:

\int_{\mathbb{R}}(g(x)-g_{N}(x))^{2}f_{X}(x)dx=\sum_{i=1}^{N}(g(x_{i})-g_{N}(x_{i}))^{2}w_{i}+\varepsilon_{N}

Donc l'erreur dans L^{2} est déterminée par l'erreur de quadrature.

Pour une variable X\sim\mathcal{N}(0,1) il existe une relation entre les pairs de \{x_{i},w_{i}\}_{i=1}^{N} et ceux donnés par la quadrature de Gauss-Hermite. En effet, la quadrature de Gauss-Hermite est basée sur la fonction de poids x\mapsto e^{-x^{2}}, pour une fonction x\mapsto\Psi(x) on approxime les intégrales de la forme \int_{-\infty}^{+\infty}e^{-x^{2}}\Psi(x)dx.

d'autre part on a \mathbb{E}(\Psi(X))=\int_{-\infty}^{+\infty}\frac{1}{\sqrt{2\pi}}e^{-\frac{x^{2}}{2}}\Psi(x)dx=\int_{-\infty}^{+\infty}\frac{1}{\sqrt{\pi}}e^{-x^{2}}\Psi(\sqrt{2}x)dx

donc la relation entre les points et les poids des deux méthodes est x_{i}^{H}=\frac{x_{i}}{\sqrt{2}}et w_{i}^{H}=w_{i}\sqrt{\pi}. L'erreur de la quadrature de Gauss-Hermite et donc de la collocation \varepsilon_{N}=\frac{N!\sqrt{\pi}}{2^{N}}\frac{\Psi^{(2N)}(\hat{\xi})}{(2N)!},\quad\Psi(x)=(g(x)-g_{N}(x))^{2}=\left(\frac{1}{N!}\frac{d^{N}g(x)}{dx^{N}}\bigg|_{x=\hat{\xi}}\prod_{i=1}^{N}(x-x_{i})\right)^{2}

Pour une fonction x\mapsto\Psi(x) assez régulière l'erreur \varepsilon_{N}converge vers 0 quand N\to\infty.

\textbf{Erreur de convergence pour les queues}

Ici on considère la différence entre Y et son approximation Y_{N} sachant que , où y^{*}détermine la queue. Pour tout i\in\llbracket1,N\rrbracket on a g_{N}(x_{i})=g(x_{i})=y_{i} et on fixe y^{*}et x^{*}=F_{X}^{-1}(F_{Y}(y^{*})). Alors dans L^{2}, on a:

\begin{array}{cc}
\mathbb{E}\left[(Y-Y_{N})^{2}|Y>y^{*}\right] & =\mathbb{E}\left[(g(X)-g_{N}(X))^{2}|X>x^{*}\right]\\
 & =\frac{1}{\mathbb{P}(X>x^{*})}\int_{-\infty}^{+\infty}(g(x)-g_{N}(x))^{2}1_{x>x*}(x)f_{X}(x)dx
\end{array}

En utilisant la quadrature

\mathbb{E}\left[(g(X)-g_{N}(X))^{2}|X>x^{*}\right]\leq\frac{1}{\mathbb{P}(X>x^{*})}\int_{-\infty}^{+\infty}(g(x)-g_{N}(x))^{2}f_{X}(x)dx=\frac{1}{\mathbb{P}(X>x^{*})}\left(\sum_{i=1}^{N}(g(x_{i})-g_{N}(x_{i}))^{2}w_{i}+\varepsilon_{N}\right)

Les deux fonctions g(x) et g_{N}(x) sont égaux dans les points de collocation donc la borne supérieure est donnée par

\mathbb{E}\left[(g(X)-g_{N}(X))^{2}|X>x^{*}\right]\leq\frac{1}{\mathbb{P}(X>x^{*})}\frac{N!\sqrt{\pi}}{2^{N}}\frac{\Psi^{(2N)}(\hat{\xi})}{(2N)!}

On peut montrer que pour x_{*}>0\mathbb{P}(X>x^{*})\geq\frac{1}{\sqrt{2\pi}}e^{-x_{*}^{2}/2}\left(\frac{1}{x^{*}}-\frac{1}{x_{*}^{3}}\right)

et pour x^{*}>1

\mathbb{E}\left[(g(X)-g_{N}(X))^{2}|X>x^{*}\right]\leq\pi\sqrt{2}e^{-x_{*}^{2}/2}\frac{x_{*}^{3}}{x_{*}^{2}-1}\frac{N!}{2^{N}(2N)!}\Psi^{(2N)}(\hat{\xi})

Donc on a \lim_{N\to\infty}\mathbb{E}\left[(g(X)-g_{N}(X))^{2}|X>x^{*}\right]=0

En utilisant l'inégalité de Chebychev on trouve

\mathbb{P}((Y-Y_{N})^{2}\geq a)\leq\frac{1}{a}\mathbb{E}((Y-Y_{N})^{2})=\frac{\varepsilon_{N}}{a}\to0

• \textbf{Elargissement de la grille de collocation}

Soit X\sim\mathcal{N}(0,1) et N le nombre de points de collocation. En calculant F_{\mathcal{N}(0,1)}(x_{i}) pour i\in\llbracket1,N\rrbracket on remarque que F_{\mathcal{N}(0,1)}(x_{1}) et F_{\mathcal{N}(0,1)}(x_{N}) sont prohce respectivement de 0 et de 1. Ceci est plus prononcé quand N est grand. Le calcul F_{Y}^{-1}(F_{\mathcal{N}(0,1)}(.)) peut révéler des instabilités numériques quand F_{\mathcal{N}(0,1)}(.)\to0 ou F_{\mathcal{N}(0,1)}(.)\to1. Pour contourner ceci, il est proposé de définir une nouvelle variable \hat{X}et calculer F_{\hat{X}}(x_{i}) à la place de F_{X}(x_{i}). On choisit la variable \hat{X}\sim\mathcal{N}(0,\sigma^{2}) tel que

F_{\mathcal{N}(0,\sigma^{2})}(x_{1})=p_{min},\quad ou\quad F_{\mathcal{N}(0,\sigma^{2})}(x_{N})=p_{max}

et donc 

\sigma=\frac{x_{1}}{F_{\mathcal{N}(0,1)}^{-1}(p_{min})}\quad ou\quad\sigma=\frac{x_{N}}{F_{\mathcal{N}(0,1)}^{-1}(p_{max})}

pour bien choisir les p_{min}et p_{max}il faut avoir une idée sur la distribution d'intérêt. Si elle a des queues lourdes il faut prendre p_{max} assez grand pour que la queue soit bien approximée par le polynôme. Une fois le \sigma déterminé l'échantillonage se fait y_{i}=F_{Y}^{-1}(F_{\mathcal{N}(0,\sigma^{2})}(x_{i})), d'où

y_{n}\approx g_{N}(\xi_{n})=\sum_{i=1}^{N}F_{Y}^{-1}\left(F_{\mathcal{N}(0,1)}\left(\frac{x_{i}}{\sigma}\right)\right)l_{i}(\xi_{n}),\quad l_{i}(\xi_{n})=\prod_{j=1,i\neq j}^{N}\frac{\sigma\xi_{n}-x_{i}}{x_{i}-x_{j}}

\textbf{Collocating local volatility model}

On considère un sous jacent S(t) et un processus (kernel process) X(t) dont on dispose des moments \mathbb{E}(X^{i}(t)), i\in\mathbb{N}. La relation entre S(t) et X(t) est donnée par 

S(t)=g(t,X(t))

où (t,x)\mapsto g(t,x) est une fonction déterministe. Le but de cette méthode est de construire la fonction g de manière à ce que les volatilités implicites générées par notre modèle soit égales à celles du marché sur un ensemble de maturités.

Le modèle est donné sous la probabilité risque neutre par

\begin{array}{c}
S(t)=g(t,X(t)),\\
dX(t)=\mu(X(t))dt+\sigma(X(t))dW^{\mathbb{Q}}(t),\quad X(t_{0})=S(t_{0})
\end{array}

Soit T_{i},i=1,\dots,M les maturités du marché. La fonction de répartition du \textit{marché}pour une maturité T_{i} est donnée par F_{\hat{S}(T_{i})}(x)=1+e^{rT_{i}}\frac{\partial C(t_{0},T_{i},K)}{\partial K}|_{K=x}, où C(t_{0},T_{i},K)=e^{-rT_{i}}\int_{K}^{\infty}(x-K)f_{\hat{S}(T_{i})}(x)dx le prix (sous probabilité risque neutre) d'une option call européene de maturité T_{i}et de strike K. On dispose aussi de la fonction de répartition du processus X(t). En utilisant la méthode de collocation décrite dans la section précédente on constuit la fonction g pour notre set de maturité et le set des points de collocation on a donc l'équation suivante avec x_{ij}:=x_{j}(T_{i}):F_{X(T_{i})}(x_{ij})=F_{\hat{S}(T_{i})}(g(T_{i},x_{ij}))=:F_{\hat{S}(T_{i})}(s_{ij}),\quad i\in\llbracket1,N\rrbracket,\quad j\in\llbracket1,M\rrbracket

les valeurs de collocation sont g(T_{i},x_{ij}):=s_{ij}=F_{\hat{S}(T_{i})}^{-1}(F_{X(T_{i})}(x_{ij})). on a donc la proprité suivante:\mathbb{E}(\halfnote\hat{S}^{k}(T_{i}))=\mathbb{E}(g^{k}(T_{i},X(T_{i}))+\varepsilon_{i,N},\quad i\in\llbracket1,M\rrbracket,\quad k\in\llbracket1,N\rrbracket

et on a montré (dans la section précédente) que l'erreur de quadrature \varepsilon_{i,N}\to0 pour la maturité T_{i} exponentiellement en N (N prend des valeurs <8).

Une fois on a la grille \{T_{i},x_{ij},s_{ij}\} on passe à l'étape suivante de la calibration en imposant la continuité à la fonction g(t,x) pour pouvoir simuler S(t) pour des maturités intérmédiaires de celles du marché (t\in]T_{i},T_{i+1}[). Il faut dans un premier temps déterminer les points de collocation x_{j}(t),\quad j\in\llbracket1,N\rrbracket,\quad t\in]T_{i},T_{i+1}[, cette procédure est décrite dans la section suivante. Pour déteminer les valeur de collocation s_{j}(t) on fait l'interpolation linéaire suivante:

\forall t\in[T_{i},T_{i+1}[,\quad s_{j}(t)=s_{ij}+(s_{i+1j}-s_{ij})\frac{t-T_{i}}{T_{i+1}-T_{i}},\forall j\in\llbracket1,N\rrbracket

Une fois les points de collocation x_{j}(t) et les valeurs de collocations s_{j}(t) déterminés, il faut déterminer une fonction continue differentiable g(t,X(t)) telle que g(t,x_{j}(t))=s_{j}(t). On utilise alors l'interpolation de Lagrange:

g(t,X(t))=\sum_{j=1}^{N}s_{j}(t)l_{j}(X(t)),\quad l_{j}(X(t))=\prod_{k=1,j\neq k}^{N}\frac{X(t)-x_{j}(t)}{x_{k}(t)-x_{j}(t)}

\textbf{Processus kernel X(t)}

Le processus X(t) est choisi librement à condition de posséder des moments. Une relation \textit{quasi linéaire}est favorable entre les variables X(T_{i}) et \hat{S}(T_{i}) pour réduire l'erreur d'approximation ou quand les densités de deux variables \textit{se ressemblent}. On peut donc considérer comme kernel process un brownien, processus qui suit la dynamique d'heston ou un Orlenstein Ulenbenk.

\textbf{Points de collocation pour une variable normale}: Soient X_{1}\sim\mathcal{N}(a_{1},b_{1}) et X_{2}\sim\mathcal{N}(a_{2},b_{2}) et leurs points de collocations respectifs x_{i}^{X_{1}} et x_{i}^{X_{2}}. Alors F_{X_{1}}(x_{i}^{X_{1}})=F_{X_{2}}(x_{i}^{X_{2}}) pour tous i\in\llbracket1,N\rrbracket et x_{i}^{X_{1}}=a_{1}+b_{1}x_{i}^{\mathcal{N}(0,1)} et x_{i}^{X_{2}}=a_{2}+b_{2}x_{i}^{\mathcal{N}(0,1)}, où x_{i}^{\mathcal{N}(0,1)}sont les points de collocation pour une variable normale standard.

En utilisant le résultat précédent, on obtient les points de collocation du processus X(t) par:

x_{i}(t)=\mathbb{E}(X(t))+\sqrt{\mathbb{V}(X(t))}x_{i}^{\mathcal{N}(0,1)},\quad i\in\llbracket1,N\rrbracket

et pour avoir les x_{i}^{\mathcal{N}(0,1)}on utilise les abscisse de Gauss-Hermite x_{i}^{H}en effet on a la relation suivante x_{i}^{\mathcal{N}(0,1)}=\sqrt{2}x_{i}^{H}.

La question qui se pose est comment choisr les paramètres du processus X(t). Considérons X_{1}(t)=X_{1}(0)+a_{1}t+b_{1}W^{\mathbb{Q}}(t) et X_{2}(t)=X_{2}(0)+a_{2}t+b_{2}W^{\mathbb{Q}}(t) (avec le même mouvement brownien) alors on a X_{2}(t)=c_{1}+c_{2}X_{1}(t) donc d'après le résultat précédent F_{X_{1}(t)}(x_{i}^{X_{1}(t)})=F_{X_{2}(t)}(x_{i}^{X_{2}(t)}) et comme la fonction g est complétement déterminée par les fonctions de répartitions alors g(X_{1}(t))=g(X_{2}(t)) p.s. Donc dans ce cas le choix des paramètres du processus X(t) n'impacte pas les résultats de la méthode de collocation. Par contre si on considère comme processus kernel un Ornstein-Uhlenbeck (OU) de dynamique dX(t)=\lambda(\theta-X(t))dt+\eta dW^{\mathbb{Q}}(t)

avec la solution:

X(t)=X_{0}e^{-\lambda t}+\theta(1-e^{-\lambda t})+\frac{\eta}{\sqrt{2\lambda}}e^{-\lambda t}W^{\mathbb{Q}}(e^{2\lambda t}-1)

ici la filtration du mouvement brownien dépend du paramètre \lambda donc si on prend deus processus OU avec \lambda_{1}\neq\lambda_{2} on aura des trajectoires g(X_{1}(t))\neq g(X_{2}(t)). 

\textbf{Least square Monte carlo}
