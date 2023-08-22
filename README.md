# particle_method
\textbf{Stochastic collocation}

Soit Y une variable aléatoire à valeurs réelles avec une fonction de répartition strictement croissante F_{Y}(y). Soit U\sim\mathcal{U}([0,1]) et u_{n}un échantillon de U. Classiquement pour générer un échantillon de Y on fait 

y_{n}=F_{Y}^{-1}(u_{n})

Or dans le cas où la fonction inverse de répartition n'a pas de forme analytique cette procédure devient couteuse.

On alors considère une autre variable X, pour laquelle F_{X}^{-1}(u_{n}) est moins couteuse que celle de Y. On sait que F_{Y}(Y)\overset{d}{=}F_{X}(X) donc y_{n}=F_{Y}^{-1}(F_{X}(\xi_{n})) où y_{n},\xi_{n} les échantillons de Y,X respectivement. Ici encore l'échantillonement de Y reste couteux. Il faut donc trouver une relation alternative pour ne pas faire l'invesion F_{Y}^{-1}pour tous l'échantillonf de X.

On cherche donc une fonction de manière à ce que g(.)=F_{Y}^{-1}(F_{X}(.)) donc à ce que Y\overset{d}{=}g(X), et telle que l'évaluation de cette fonction n'est pas couteuse.

Dans la méthode de collocation stochastique on approxime Y par une foncition g de X en terme d'expansion de Lagrange l_{i}(\xi_{n}):y_{n}\approx g_{N}(\xi_{n})=\sum_{i=1}^{N}y_{i}l_{i}(\xi_{n}),\quad l_{i}(\xi_{n})=\prod_{j=1,i\neq j}^{N}\frac{\xi_{n}-x_{j}}{x_{i}-x_{j}}

où \xi_{n}est un échantillon de X et x_{i},x_{j}sont des points de collocation, et y_{i}=F_{Y}^{-1}(F_{X}(x_{i})). \textbf{l}(x)=(l_{1}(x),\dots,l_{N}(x))^{T} est la base de Lagrange., telle que l_{i}(x_{j})=\delta_{ij}. Donc une fois les N points de collocation déterminés x_{i} et les N inversions F_{Y}^{-1} faites, on peut simuler nimporte quel nombre d'échantillons de Y et ceci par l'évaluation du polynome g_{N}(.). On parle ici de \textit{Stochastic Collocation Monte Carlo sampler}.

\textbf{Points de collocation}

On dit qu'une sequence de polynomes orthogonaux \{p_{i}\}_{i=0}^{N} avec degré deg(p_{i})=i est orthogonale en L^{2}par rapport à la densité de f_{X}(X) de X, si

\mathbb{E}(p_{i}(X)p_{j}(X))=\delta_{ij}\mathbb{E}(p_{i}^{2}(X)),\quad i,j=0,\dots,N

Et pour toutes densité f_{X}(.), il existe une sequence de polynomes orthogonaux p_{i}(x) unique avec comme degré deg(p_{i}(x))=i, cette sequence se construit par 

p_{i+1}(x)=(x-\alpha_{i})p_{i}(x)-\beta_{i}p_{i-1}(x),\quad i\in\llbracket0,N-1\rrbracket

où p_{-1}(x)=0 et p_{0}(x)=1 et pour tous i\in\llbracket0,N-1\rrbracket \alpha_{i}=\frac{\mathbb{E}(Xp_{i}^{2}(X))}{\mathbb{E}(p_{i}^{2}(X))} et \beta_{i}=\frac{\mathbb{E}(p_{i}^{2}(X))}{\mathbb{E}(p_{i-1}^{2}(X))} et \beta_{0}=0

On construit la matrice de Gram M=\{\mu_{ij}\}_{i,j=0}^{N} en considérant le monomial m_{i}(X)=X^{i} par \mu_{ij}=\mathbb{E}(m_{i}(X)m_{j}(X))=\mathbb{E}(X^{i+j}). La matrice M est définie positive elle s'écrit alors comme M=R^{T}R où R=\{r_{ij}\}_{i,j=0}^{N} est une matrice triangulaire inférieure. On a 

\alpha_{j}=\frac{r_{j,j+1}}{r_{j,j}}-\frac{r_{j-1,j}}{r_{j-1,j-1}},\quad\beta_{j}=\left(\frac{r_{j+1,j+1}}{r_{j,j}}\right)^{2},\quad j=1,\dots,N-1

où r_{0,0}=1 et r_{0,1}=0

Les zeros x_{i},i\in\llbracket1,N\rrbracket du polynome orthogonal p_{N}(X) sont les valeurs propres de la matrice symétrique suivante

\hat{J}:=\left(\begin{array}{ccccc}
\alpha_{1} & \sqrt{\beta_{1}} & 0 & 0 & 0\\
\sqrt{\beta_{1}} & \alpha_{2} & \sqrt{\beta_{2}} & 0 & 0\\
 & \ddots & \ddots & \ddots\\
0 & 0 & \sqrt{\beta_{N-2}} & \alpha_{N-1} & \sqrt{\beta_{N-1}}\\
0 & 0 & 0 & \sqrt{\beta_{N-1}} & \alpha_{N}
\end{array}\right)

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
