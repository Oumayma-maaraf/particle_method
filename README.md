\textbf{Problèmatique}

On considère un portefeuille de transactions entre une banque et une contrepartie collatéralisée générant entre deux temps.

On considère un portefeuille de transactions entre une banque et une contrepartie collatéralisée générant \Pi(t,T) entre t et la maturité T. La valeur de ce portefeuille en t est V_{t}=\mathbb{E}(\Pi(t,T)|\mathcal{F}_{t})

On note \tau_{c}le temps de défault de la contrepartiertie et R_{c}le taux de recouvrement dans le cas du défault. La CVA est l'espérence des pertes futures dues aux défauts conditionnel aux états de crédit. On peut la voir comme la différence entre le prix du produit avec et sans prendre en considération le défaut de la contrepartie. En faisant l'hypothèse de l'indépendance entre le temps de défaut et la valeur du portefeuille, on trouve CVA(t)=-(1-R_{c})\int_{t}^{T}\mathbb{E}\left[B(t,T)(V_{s}-C_{s})_{+}|\mathcal{F}_{t}\right]d\mathbb{Q}(\tau_{c}<T)

Donc pour t=0 on a 

CVA(0)=-(1-R_{c})\int_{0}^{T}\mathbb{E}\left[B(0,T)(V_{s}-C_{s})_{+}\right]d\mathbb{Q}(\tau_{c}<T)

Où B(0,t)=\exp\left(-\int_{0}^{t}r_{s}ds\right) est le facteur d'actualisation et C_{t} le collatéral à la date t en effet on a C_{t}=f(V_{t}) avec f une fonction déterministe.

Considérant t_{0}=t<t_{1}<\dots<t_{N}=T, en pratique la CVA est calculée par 

CVA(t)=-(1-R_{c})\sum_{k=1}^{N}\mathbb{E}\left[B(t,t_{k})(V_{t_{k}}-C_{t_{k}})_{+}\bigg|\mathcal{F}_{t}\right]\mathbb{Q}(t_{k-1}<\tau_{c}<t_{k})

Dans le but d'augmenter la précision du calcul du CVA(0) tout en gardant le nombre de simulation des sous jacents constant (N), on peut soit s'intéresser à \mathbb{E}\left[B(0,T)(V_{t}-C_{t})_{+}\right] et faire des techniques de réduction de variance, soit s'intéresser à V_{t}=\mathbb{E}(\Pi(t,T)|\mathcal{F}_{t}) et essayer d'augmenter la précision de son calcul. Dans le dernier cas, il peut choisir pour notre actif sous jacent une dynamique qui modélise bien son cours dans le marché et les prix des options qui ont dépendent. Or qu'en xVa on travaille sur plusieurs asset classe en même temps avec des structures de corrélations assez complexes et on ne se permet pas de faire pour chaque asset classe un modèle sophistiqué. Du coup ce qu'on propose dans ce travail est une approche modèle free dans laquelle (ou dépend d'un modèle simple) afin de reproduire les résultats pour plusieurs asset classes. 

\textbf{Modèle de marché}

On considère que le marché est donné par un modèle qu'on spécifiera. Soit (\Omega,(\mathcal{F}_{t})_{0\leq t},\mathbb{P}) un espace de probabilité filtré, soit W_{t} un mouvement brownien sous \mathbb{P}. L'actif sans risque suit la dynamique suivantedB_{t}=rB_{t}dt

On considère un modèle à volatilité locale donné par le modèle CEV (Constant elasticity of variance). Dans ce modèle la volatilité est une fonction puissance du niveau de sous jacent:\frac{dS_{t}}{S_{t}}=\mu dt+\sigma_{0}S_{t}^{\alpha-1}dW_{t},

Le modèle de Black-Scholes et le modèle gaussien sont des cas limites de cette équation vec \alpha=1 et \alpha\to0 respectivement. On suppose que le prix forward du sous jacent F_{t}=e^{r(T-t)}S_{t}suit le modèle CEV

dF_{t}=\sigma_{0}F_{t}^{\alpha}dW_{t},\quad0<\alpha\leq1

La valeur 0 est une barrière absorbante, si F_{t}=0 pour un t, alors F_{s}=0 pour tout s\geq t.

\textbf{Propriété de martingale}

On montre que pour 0<\alpha\leq1, F_{t} définit une vraie martiangle de carré intégrable sur [0,T] pour tout T<\infty. Soit \tau_{n}=\inf\{t:F_{t}\geq n\}, F_{T\wedge\tau_{n}}est alors de carré intégrable, et on a pour 0<\alpha\leq1\mathbb{E}[F_{\tau_{n}\wedge T}^{2}]=\sigma_{0}\mathbb{E}\left[\int_{0}^{\tau_{n}\wedge T}F_{t}^{2\alpha}dt\right]\leq\sigma_{0}^{2}\mathbb{E}\left[\int_{0}^{\tau_{n}\wedge T}(1+F_{t}^{2})dt\right]\leq\sigma_{0}^{2}\mathbb{E}\left[\int_{0}^{T}(1+F_{t\wedge\tau_{n}}^{2})dt\right]

Par le lemme de Gronwall on a:

\mathbb{E}\left[F_{\tau_{n}\wedge T}^{2}\right]\leq\sigma_{0}^{2}Te^{\sigma_{0}^{2}T}

d'où par convergence monotone\mathbb{E}\left\{ \sigma_{0}^{2}\int_{0}^{T}F_{t}^{2\alpha}dt\right\} <\infty

\textbf{Volatilité implicite}

La forme de la volatilité implicite du modèle CEV est connue grâce à l'approximation asymptotique de Hagan et Woodward:

\sigma^{imp}(K,T)=\frac{\sigma_{0}}{F_{m}^{1-\alpha}}\left\{ 1+\frac{(1-\alpha)(2+\alpha)}{24}\left(\frac{F_{0}-K}{F_{m}}\right)^{2}+\frac{(1-\alpha)^{2}}{24}\frac{\sigma_{0}^{2}T}{F_{m}^{2-2\alpha}}+\dots\right\} ,\quad F_{m}=\frac{1}{2}(F_{0}+K)

au premier ordre, on a donc \sigma^{imp}(K,T)\approx\frac{\sigma_{0}}{F_{m}^{1-\alpha}}: la volatilité implicite a la même forme que la volatilité locale mais avec une pente à la monnaie plus petite.

\textbf{Aspect pratique}

Pour utiliser ce modèle en tant que marché, on diffuse notre sous jacent S_{t} on génère une surface de prix pour plusieurs maturités et strikes après on implicite les volatilités on se retouve donc avec une surface de volatilité.

On introduit maintenant la première méthode \textit{empirique}dans laquelle on essaye de reconstruire la distribution de F_{S_{T}}ici je vais insérer des graphs pour illustrer la densité et comment on la reconstruit faire les traits de K_i 

\textbf{Déformation d'échantillon}

Soit T une maturité, on suit la procédure suivante 

• On cherche 0=K_{0}<K_{1}<\dots<K_{n}<K_{n+1}=+\infty tels que \mathbb{Q}(S_{T}\in[K_{i},K_{i+1}[)=\frac{1}{n+1}, càd de résoudre le système d'équations suivant \begin{cases}
F_{S_{T}}(K_{1})-F_{S_{T}}(K_{0}) & =\frac{1}{n+1}\\
\vdots & \vdots\\
1-F_{S_{T}}(K_{n}) & =\frac{1}{n+1}
\end{cases}

• Simuler \left(S_{T}^{i}\right)_{0\leq i\leq N}et les ordonner \left(\tilde{S}_{T}^{i}\right)_{0\leq i\leq N}

• On construit n bucket (B_{i}) tels que chaque bucket contient \left\lfloor \frac{N}{n}\right\rfloor simulations ordonnées. 

• On définie \forall i\in\llbracket0,n\rrbracket,\quad\mathcal{I}_{i}=\{\mbox{indexes of}\quad\tilde{S}_{T}\in B_{i}\}et Card_{B}=\#B. On assure que \forall i\in\llbracket0,n\rrbracket, j\in\mathcal{I}_{i} \tilde{S}_{T}^{j}\in[K_{i},K_{i+1}]. Pour ce faire on fait la transformation suivante sur tous les éléments de chaque bucket i\in\llbracket0,n-1\rrbracket\begin{cases}
\tilde{S}_{T}^{0*} & =K_{i}\\
\tilde{S}_{T}^{1*} & =\tilde{S}_{T}^{0*}+\frac{\tilde{S}_{T}^{n*}-\tilde{S}_{T}^{0*}}{\tilde{S}_{T}^{n}-\tilde{S}_{T}^{0}}(\tilde{S}_{T}^{1}-\tilde{S}_{T}^{0})\\
\vdots & \vdots\\
\tilde{S}_{T}^{n*} & =K_{i+1}
\end{cases}

Pour le dernier bucket on fait la transformation suivante \forall j\in\mathcal{I}_{n},\quad\tilde{S}_{t}^{j*}=\max(\tilde{S}_{t}^{j},K_{n})

• On résoud le système d'équation suivant (on cherche \forall i\in\llbracket n,0\rrbracket \alpha_{i}), on note Num le numéraire \frac{1}{N}\sum_{l\in\bigcup_{j=1}^{n}\mathcal{I}_{j}}\frac{(\hat{S}_{T}^{l*}-K_{i})_{+}}{Num^{l}}=\mathcal{P}^{theo}(C(T,K_{i}))

On prend K_{n+1}=\max(N\mathcal{P}^{theo}(C(T,K_{i}))+K_{n},\max(\tilde{S}_{T}^{n*})). Si l\in\mathcal{I}_{i}\hat{S}_{T}^{l*}=\begin{cases}
\alpha_{i}\tilde{S}_{T}^{l*}+(1-\alpha_{i})K_{i} & si\quad\frac{1}{N}\left(\sum_{l\in\mathcal{I}_{i}}\frac{(\tilde{S}_{T}^{l*}-K_{i})_{+}}{Num^{l}}+\sum_{l\in\bigcup_{j=i+1}^{n}\mathcal{I}_{j}}\frac{(\hat{S}_{T}^{l*}-K_{i})+}{Num^{l}}\right)\geq\mathcal{P}^{theo}(C(T,K_{i}))\\
\alpha_{i}\tilde{S}_{T}^{l*}+(1-\alpha_{i})K_{i+1} & si\quad\frac{1}{N}\left(\sum_{l\in\mathcal{I}_{i}}\frac{(\tilde{S}_{T}^{l*}-K_{i})_{+}}{Num^{l}}+\sum_{l\in\bigcup_{j=i+1}^{n}\mathcal{I}_{j}}\frac{(\hat{S}_{T}^{l*}-K_{i})+}{Num^{l}}\right)\leq\mathcal{P}^{theo}(C(T,K_{i}))
\end{cases}

On peut résoudre de manière itérative l'équation précédente où on note x\in\{i,i+1\}:

how to pass time exactly how to pass two hours how i wil show do something anything o r nothing just pretend to be super 

\sum_{l\in\bigcup_{j=i}^{n}}\frac{(\hat{S}_{T}^{l*}-K_{i})}{Num^{l}}=N\mathcal{P}^{theo}(C(T,K_{i}))

\sum_{l\in\mathcal{I}_{i}}\frac{\hat{S}_{T}^{l*}}{Num^{l}}-K_{i}\sum_{l\in\mathcal{I}_{i}}\frac{1}{Num^{l}}-K_{i}\sum_{k=i+1}^{n}\sum_{l\in\mathcal{I}_{k}}\frac{1}{Num^{l}}+\sum_{k=i+1}^{n}\sum_{l\in\mathcal{\mathcal{I}}_{k}}\frac{\hat{S}_{T}^{l*}}{Num^{l}}=N\mathcal{P}^{theo}(C(T,K_{i}))

\textbf{Stochastic collocation}

Soit Y une variable aléatoire à valeurs réelles et F_{Y}(y) sa fonction de répartition qui est strictement croissante . Soit U\sim\mathcal{U}([0,1]) et u_{n} un échantillon de U. Classiquement pour générer un échantillon de Y on fait 

y_{n}=F_{Y}^{-1}(u_{n})

Or dans le cas où la fonction inverse de répartition n'a pas de forme analytique cette procédure devient couteuse, car il faut faire autant d'inversion que d'échantillon.

On considère alors une autre variable X, pour laquelle F_{X}^{-1}(.) est moins couteuse que celle de Y. On sait que F_{Y}(Y)\overset{d}{=}F_{X}(X) donc y_{n}=F_{Y}^{-1}(F_{X}(\xi_{n})) où y_{n},\xi_{n} les échantillons de Y,X respectivement. Ici encore l'échantillonement de Y reste couteux, car on fait autant d'inversion que d'échantillon. Il faut donc trouver une relation alternative pour ne pas faire l'invesion F_{Y}^{-1}pour tout l'échantillon de X.

On cherche donc une fonction g de manière à ce que g(.)=F_{Y}^{-1}(F_{X}(.)) donc à ce que Y\overset{d}{=}g(X), et telle que l'évaluation de cette fonction n'est pas couteuse.

Dans la méthode de collocation stochastique on approxime Y par une fonction g de X en terme d'expansion de Lagrange l_{i}(\xi_{n}):y_{n}\approx g_{N}(\xi_{n})=\sum_{i=1}^{N}y_{i}l_{i}(\xi_{n}),\quad l_{i}(\xi_{n})=\prod_{j=1,i\neq j}^{N}\frac{\xi_{n}-x_{j}}{x_{i}-x_{j}}

où \xi_{n}est un échantillon de X et x_{i},x_{j}sont des points de collocation (N est généralement <8), et y_{i}=F_{Y}^{-1}(F_{X}(x_{i})). \textbf{l}(x)=(l_{1}(x),\dots,l_{N}(x))^{T} est la base de Lagrange., telle que l_{i}(x_{j})=\delta_{ij}. Donc une fois les N points de collocation déterminés x_{i} et les N inversions F_{Y}^{-1} faites, on peut simuler nimporte quel nombre d'échantillons de Y et ceci par l'évaluation du polynome g_{N}(.). On parle ici de \textit{Stochastic Collocation Monte Carlo sampler}.

• \textbf{Points de collocation}

Les points de collocation optimaux sont choisis pour être les points de la quadrature de Gauss et sont définis comme les zéros du polynôme orthogonal correspondant. On sait que g_{N}(x)=\sum_{i=1}^{N}y_{i}l_{i}(x),\quad l_{i}(x)=\prod_{j=1,i\neq j}^{N}\frac{x-x_{j}}{x_{i}-x_{j}}

Pour éviter de faire O(N^{2}) opération pour chaque nouvelle valeurs de x, on considère un poids \lambda_{i} défini par 

\lambda_{i}=\frac{1}{\prod_{j=1,j\neq i}^{N}(x_{i}-x_{j})},\quad i\in\llbracket1,N\rrbracket

et l(x)=(x-x_{1})\dots(x-x_{N}), donc l_{i}(x) peut s'écrire comme

l_{i}(x)=l(x)\frac{\lambda_{i}}{x-x_{i}}

donc g_{N}(x)=\sum_{i=1}^{N}\frac{y_{i}\lambda_{i}}{(x-x_{i})}l(x)

On dit qu'une sequence de polynomes orthogonaux \{p_{i}\}_{i=0}^{N} avec degré deg(p_{i})=i est orthogonale en L^{2}par rapport à la densité de f_{X}(X) de X, si

\mathbb{E}\left[p_{i}(X)p_{j}(X)\right]=\delta_{ij}\mathbb{E}\left[p_{i}^{2}(X)\right],\quad i,j=0,\dots,N

Et pour toute densité f_{X}(.), il existe une sequence de polynomes orthogonaux p_{i}(x) unique avec comme degré deg(p_{i}(x))=i, cette séquence se construit par 

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

Les points de collocation pour une variable X\sim\mathcal{N}(0,1) sont donnés par 

\begin{array}{cccccccccccc}
x_{i} & N=2 & N=3 & N=4 & N=5 & N=6 & N=7 & N=8 & N=9 & N=10 & N=11\\
x_{1} & -1 & -1.7321 & -2.3344 & -2.8570 & -3.3243 & -3.7504 & -4.1445 & -4.5127 & -4.8595 & -5.1880\\
x_{2} & 1 & 0.0 & -0.7420 & -1.3556 & -1.8892 & -2.3668 & -2.8025 & -3.2054 & -3.5818 & -3.9362\\
x_{3} &  & 1.7321 & 0.7420 & 0.0 & -0.6167 & -1.1544 & -1.6365 & -2.0768 & -2.4843 & -2.8651\\
x_{4} &  &  & 2.3344 & 1.3556 & 0.6167 & 0.0 & -0.5391 & -1.0233 & -1.4660 & -1.8760\\
x_{5} &  &  &  & 2.8570 & 1.8892 & 1.1544 & 0.5391 & 0.0 & -0.4849 & -0.9289\\
x_{6} &  &  &  &  & 3.3243 & 2.3668 & 1.6365 & 1.0233 & 0.4849 & 0.0\\
x_{7} &  &  &  &  &  & 3.7504 & 2.8025 & 2.0768 & 1.4660 & 0.9289\\
x_{8} &  &  &  &  &  &  & 4.1445 & 3.2054 & 2.4843 & 1.8760\\
x_{9} &  &  &  &  &  &  &  & 4.5127 & 3.5818 & 2.8651\\
x_{10} &  &  &  &  &  &  &  &  & 4.8595 & 3.9362\\
x_{11} &  &  &  &  &  &  &  &  &  & 3.9362
\end{array}

\textbf{Monotonie}

Une fois les points de collocations déterminés x_{i} et les inversions correspondantes faites y_{i}=F_{Y}^{-1}(F_{X}(x_{i})) on doit construire une fonction d'approximation g_{N}(x) qui est \textit{idéalement}monotone différentiable et qui vérifie y_{i}=g_{N}(x_{i}). en effet en choisissant g_{N}(x) comme un polynôme de Lagrange en ne garantie pas la monotonie. Néanmoins la convergence du \textit{SCMC sampler} ne dépend pas sur la monotonie de g_{N}(x).

• \textbf{Analyse d'erreur}

On s'intéresse dans cette section à l'erreur générée par \textit{Stochastic Collocation Monte Carlo sampler}.

On se met dans un premier temps dans un cas où la méthode de collocation donne des résultats exacts. Soient Y\sim\mathcal{N}(\mu_{Y},\sigma_{Y}^{2}) et X\sim\mathcal{N}(\mu_{X},\sigma_{X}^{2}) deux variables aléatoires alors g_{N}(X)\overset{d}{=}Y pour N=2. En effet, soietn x_{1}et x_{2} deux points de collocation alorsg_{2}(X)=y_{1}\frac{X-x_{2}}{x_{1}-x_{2}}+y_{2}\frac{X-x_{1}}{x_{2}-x_{1}}

On a F_{\mathcal{N}(0,1)}\left(\frac{y_{i}-\mu_{Y}}{\sigma_{Y}}\right)=F_{\mathcal{N}(0,1)}\left(\frac{x_{i}-\mu_{X}}{\sigma_{X}}\right), alors y_{i}=\frac{x_{i}-\mu_{X}}{\sigma_{X}}\sigma_{Y}+\mu_{Y}. On a donc \mathbb{E}\left(g_{2}(X)\right)=\mu_{Y} et \mathbb{V}\left(g_{2}(X)\right)=\sigma_{Y}^{2} et comme g_{2}(X) suit une loi normale alors Y\overset{d}{=}g_{2}(X).

Dans un cas général, pour mesurer l'erreur on peut soit considérer la différence entre g(X) et g_{N}(X) soit l'erreur associée à l'approximation de la fonction de répartition.

La première erreur est liée à l'interpolation de Lagrange, en effet la relation entre Y et X est Y=g(X) qu'on approxime par un polynome de Lagrange Y\approx g_{N}(X) pour N points de collocation. Cette erreur est donc bien connue 

e_{X}(\xi_{n})=|g(\xi_{n})-g_{N}(\xi_{n})|=\bigg|\frac{1}{N!}\frac{d^{N}g(x)}{dx^{N}}\bigg|_{x=\hat{\xi}}\prod_{i=1}^{N}(\xi_{n}-x_{i})\bigg|

avec x_{i}est un point de collocation, \hat{\xi}\in[\min(x),\max(x)] et x=(x_{1},\dots,x_{N})^{T}, on peut borner cette erreur en prenant \hat{\xi} l'abscisse du maximum de \bigg|\frac{d^{N}g(x)}{dx^{N}}\bigg|. En utilisant\xi_{n}=F_{X}^{-1}(u_{n}), on trouve

e_{U}(u_{n})=\bigg|g\left(F_{X}^{-1}(u_{n})\right)-g_{N}\left(F_{X}^{-1}(u_{n})\right)\bigg|=\bigg|\frac{1}{N!}\frac{d^{N}g(x)}{dx^{N}}\bigg|_{x=\hat{\xi}}\prod_{i=1}^{N}\left(F_{X}^{-1}(u_{n})-x_{i}\right)\bigg|.

\textbf{Erreur de convergence en}L^{2}

On a Y=g(X)\approx Y_{N}\equiv g_{N}(X), où g(x)=F_{Y}^{-1}(F_{X}(x)) donc 

\mathbb{E}\left[(Y-Y_{N})^{2}\right]=\mathbb{E}\left[(g(X)-g_{N}(X))^{2}\right]=\int_{\mathbb{R}}(g(x)-g_{N}(x))^{2}f_{X}(x)dx

Les points de collocations x_{i} et les poinds w_{i}sont déterminés par le théorème A.2. Comme g(x_{i})=g_{N}(x_{i}),pour i\in\llbracket1,N\rrbracket l'erreur est:

\int_{\mathbb{R}}(g(x)-g_{N}(x))^{2}f_{X}(x)dx=\sum_{i=1}^{N}(g(x_{i})-g_{N}(x_{i}))^{2}w_{i}+\varepsilon_{N}=\varepsilon_{N}

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

Une fois on a la grille \{T_{i},x_{ij},s_{ij}\} on passe à l'étape suivante de la calibration en imposant la continuité à la fonction g(t,x) pour pouvoir simuler S(t) pour des maturités autres que celles du marché (t\in]T_{i},T_{i+1}[). Il faut dans un premier temps déterminer les points de collocation x_{j}(t),\quad j\in\llbracket1,N\rrbracket,\quad t\in]T_{i},T_{i+1}[, cette procédure est décrite dans la section suivante. Pour déteminer les valeur de collocation s_{j}(t) on fait l'interpolation linéaire suivante:

\forall t\in[T_{i},T_{i+1}[,\quad s_{j}(t)=s_{ij}+(s_{i+1j}-s_{ij})\frac{t-T_{i}}{T_{i+1}-T_{i}},\forall j\in\llbracket1,N\rrbracket

Une fois que les points de collocation x_{j}(t) et les valeurs de collocations s_{j}(t) sont déterminés, il faut déterminer une fonction continue differentiable g(t,X(t)) telle que g(t,x_{j}(t))=s_{j}(t). On utilise alors l'interpolation de Lagrange:

g(t,X(t))=\sum_{j=1}^{N}s_{j}(t)l_{j}(X(t)),\quad l_{j}(X(t))=\prod_{k=1,j\neq k}^{N}\frac{X(t)-x_{j}(t)}{x_{k}(t)-x_{j}(t)}

\textbf{Processus kernel X(t)}

Le processus X(t) est choisi librement à condition de posséder des moments. Une relation \textit{quasi linéaire}est favorable entre les variables X(T_{i}) et \hat{S}(T_{i}) pour réduire l'erreur d'approximation. On peut donc considérer comme kernel process un brownien, processus qui suit la dynamique d'heston ou un Orlenstein Ulenbenk.

\textbf{Points de collocation pour une variable normale}: Soient X_{1}\sim\mathcal{N}(a_{1},b_{1}) et X_{2}\sim\mathcal{N}(a_{2},b_{2}) et leurs points de collocations respectifs x_{i}^{X_{1}} et x_{i}^{X_{2}}. Alors F_{X_{1}}(x_{i}^{X_{1}})=F_{X_{2}}(x_{i}^{X_{2}}) pour tous i\in\llbracket1,N\rrbracket et x_{i}^{X_{1}}=a_{1}+b_{1}x_{i}^{\mathcal{N}(0,1)} et x_{i}^{X_{2}}=a_{2}+b_{2}x_{i}^{\mathcal{N}(0,1)}, où x_{i}^{\mathcal{N}(0,1)}sont les points de collocation pour une variable normale standard.

En utilisant le résultat précédent, on obtient les points de collocation du processus X(t) par:

x_{i}(t)=\mathbb{E}(X(t))+\sqrt{\mathbb{V}(X(t))}x_{i}^{\mathcal{N}(0,1)},\quad i\in\llbracket1,N\rrbracket

et pour avoir les x_{i}^{\mathcal{N}(0,1)}on utilise les abscisses de Gauss-Hermite x_{i}^{H}en effet on a la relation suivante x_{i}^{\mathcal{N}(0,1)}=\sqrt{2}x_{i}^{H}.

La question qui se pose est comment choisir les paramètres du processus X(t). Considérons X_{1}(t)=X_{1}(0)+a_{1}t+b_{1}W^{\mathbb{Q}}(t) et X_{2}(t)=X_{2}(0)+a_{2}t+b_{2}W^{\mathbb{Q}}(t) (avec le même mouvement brownien) alors on a X_{2}(t)=c_{1}+c_{2}X_{1}(t) donc d'après le résultat précédent F_{X_{1}(t)}(x_{i}^{X_{1}(t)})=F_{X_{2}(t)}(x_{i}^{X_{2}(t)}) et comme la fonction g est complétement déterminée par les fonctions de répartitions alors g(X_{1}(t))=g(X_{2}(t)) p.s. Donc dans ce cas le choix des paramètres du processus X(t) n'impacte pas les résultats de la méthode de collocation. Par contre si on considère comme processus kernel un Ornstein-Uhlenbeck (OU) de dynamique dX(t)=\lambda(\theta-X(t))dt+\eta dW^{\mathbb{Q}}(t)

avec la solution:

X(t)=X_{0}e^{-\lambda t}+\theta(1-e^{-\lambda t})+\frac{\eta}{\sqrt{2\lambda}}e^{-\lambda t}W^{\mathbb{Q}}(e^{2\lambda t}-1)

ici la filtration du mouvement brownien dépend du paramètre \lambda donc si on prend deus processus OU avec \lambda_{1}\neq\lambda_{2} on aura des trajectoires g(X_{1}(t))\neq g(X_{2}(t)). 

\textbf{Least square Monte carlo}

Soient T>0 un horizon de temps fixé et \mathbb{Q}l'unique probabilité risque neutre. On suppose que pour tout w\in\Omega, la trajectoire \tilde{w}:t\mapsto S_{t}(w) à valeurs dans \mathbb{R}^{d} appartient à un certain espace fonctionnelle que l'on note \mathbb{D}_{d}([0,T]):=\mathbb{D}_{d}([0,T],\mathbb{R}^{d}).

Soit F_{T}:\mathbb{D}_{d}([0,T))\mapsto\mathbb{R}^{l_{T}}une fonctionnelle mesurable. Soit g:\mathbb{R}^{l_{T}}\mapsto\mathbb{R}une fonction borélienne.

On note X le payoff stochastique d'un actif contingent. X est \mathcal{F}_{T} mesurable. Le prix de l'actif contingent est défini, à tout instant t\in[0,T] par:

\mathbb{E}_{\mathbb{Q}}\left(X\exp\left(-\int_{t}^{T}r_{s}ds\right)\bigg|\mathcal{F}_{t}\right)

On suppose que (r_{s})_{s}est constant à tout instant et on s'intéresse particulièrement à \mathbb{E}_{\mathbb{Q}}(X|\mathcal{F}_{t})

On suppose que X\in L^{2}(\mathbb{Q}), on a \mathbb{E}_{\mathbb{Q}}(X|\mathcal{F}_{t})=g(F_{t}(S))

donc g\in L^{2}\left(\mathbb{R}^{l_{t}},\mathcal{B}(\mathbb{R}^{l_{t}}),\mathbb{Q}_{F_{t}(Z)}\right)et on note L_{t}^{2}:=L^{2}\left(\mathbb{R}^{l_{t}},\mathcal{B}(\mathbb{R}^{l_{t}}),\mathbb{Q}_{F_{t}(Z)}\right)

On sait que L_{t}^{2}muni du produit scalaire:

\langle f,h\rangle_{t}=\int_{\mathbb{R}^{l_{t}}}f(x)g(x)d\mathbb{Q}_{F_{t}(Z)}(x)=\mathbb{E_{Q}}\left[f(F_{t}(S))h(F_{t}(S))\right]

est un espace d'hilbert séparable. Elle admet donc une base hilbertienne dénombrable:\left(\exists(\nu_{k})_{k\geq0}\in(L_{t}^{2})^{\mathbb{N}}orthonormale\right)\quad\left(\exists(\beta_{k})_{k\geq0}\in\mathbb{R}^{\mathbb{N}}\right):\quad g=\sum_{k=0}^{+\infty}\beta_{k}\nu_{k}

tel que

(\forall k\in\mathbb{N}):\quad\beta_{k}=\langle g,\nu_{k}\rangle=\mathbb{E_{Q}}\left[g\left(F_{t}(S)\right)\nu_{k}\left(F_{t}(S)\right)\right]=\mathbb{E_{Q}}\left[\mathbb{E_{Q}}(X|\mathcal{F}_{t})\nu_{k}\left(F_{t}(S)\right)\right]=\mathbb{E_{Q}}\left[\mathbb{E_{Q}}(X\nu_{k}\left(F_{t}(S)\right)|\mathcal{F}_{t}\right]

ainsi:

(\forall k\in\mathbb{N}):\quad\beta_{k}=\mathbb{E_{Q}}(X\nu_{k}(F_{t}(S))

On définit l'erreur de projection par p_{t}(F_{T}(S)):=X-g(F_{T}(S))

La méthode LSMC essaye d'estimer la fonction g à partir de son éciture dans la base hilberienne. En utilisant des données simulées sous \mathbb{Q}. Toutefois, en pratique, on ne peut pas simuler une infinité de données. Sur ce on essaye de travailler sur une sous-famille finie de (\nu_{k})_{k}que l'on note (\nu_{k})_{1\leq k\leq K}.

Ainsi on approxime g par g_{K}définie parg_{K}:=\sum_{k=1}^{K}\beta_{k}\nu_{k}

l'erreur d'approximation est définie par:

a_{K}:=g-g_{K}

On a a_{K}\underset{K\to+\infty}{\to}0 et par orthogonalité de la base on a\langle g_{K},a_{K}\rangle=\mathbb{E_{Q}}[g_{K}(F_{t}(S))a_{K}(F_{t}(S))]=0

On peut maintenant écrire la régression suiavnte:X=g_{K}(F_{t}(S))+a_{K}(F_{t}(S))+p_{t}(F_{T}(S))

Maintenant, étant donné un échantillon simulé de taille N: \left(\left(x^{j},A_{t}(S^{j})\right)\right)_{1\leq j\leq N} Il est naturel d'estimer g_{K}par la projection:

\hat{g}_{K}:=\arg\underset{g\in\mathcal{H}_{K}}{\min}\frac{1}{N}\sum_{j=1}^{N}\left(x^{j}-g\left(A_{t}(S^{j})\right)\right)^{2}

tel que

\mathcal{H}_{K}:=\left\{ g:\mathbb{R}^{l_{t}}\mapsto\mathbb{R}:\left(\exists(\beta_{1},\dots,\beta_{K})\in\mathbb{R}^{K}\right),g=\sum_{k=1}^{N}\beta_{k}\nu_{k}\right\} 

\textbf{Les polynômes orthogonaux de Laguerre}:

Soit n\in\mathbb{N}, on définit les fonctions \psi_{n} et L_{n}par

\left(\forall t\in\mathbb{R}\right):\quad\psi_{n}(t)=e^{-t}\frac{t^{n}}{n!};\quad L_{n}(t)=\psi_{n}^{(n)}(t)e^{t}

La famille \left(L_{n}\right)_{n}est une famille de polynôme de deg n appelé: polynôme de Laguerre. On peut démontrer que:\left(\forall n\in\mathbb{N}\right):\quad L_{n}:=\sum_{k=0}^{n}\mbox{C}_{n}^{k}\frac{(-1)^{k}}{k!}X^{k}

Utilisant ces polynômes, on désire construire une base orthonormale de l'espace L^{2}(\mathbb{R}_{+}). On définit alors la famille des fonctions (\phi_{n})_{n}, par:

\left(\forall t\in\mathbb{R}_{+}\right):\quad\phi_{n}(t):=e^{-\frac{t}{2}}L_{n}(t)

On a bien \left(\phi_{n}\right)_{n}\in L^{2}(\mathbb{R}_{+})^{\mathbb{N}}et 

\left(\forall n,p\in\mathbb{N}\right):\quad\langle\phi_{n},\phi_{p}\rangle=\int_{0}^{+\infty}L_{n}(t)L_{p}(t)e^{-t}dt=\delta_{n,p}

La dernière égalité est obtenur par l'othonormalité de la famille des polynômes de Laguerre par rapport aux produit scalaire défini sur \mathbb{R}[X], par:

\langle P,Q\rangle:=\int_{0}^{+\infty}P(t)Q(t)e^{-t}dt

On se place maintenant dans le cadre du modèle de Black-Scholes-Merton classique (d=1) Considérons le problème du pricing d'une option asiatique arithmétique et discrète:

(t_{j})_{0\leq j\leq K}une subdivision ordonnée de [0,T], tel que t_{0}=0 et t_{K}=TX:=\left(\left(1+K\right)^{-1}\sum_{j=0}^{N}S_{t_{j}}-K\right)_{+}

On a (\forall\tilde{w}\in D_{1}):\quad A_{t_{j}}(\tilde{w}):=(j+1)^{-1}\sum_{i=0}^{j}\tilde{w}(t_{i})

Maintenant on observe que l_{t_{j}}=1 et g_{t_{j}}:\mathbb{R}_{+}\to\mathbb{R}_{+}g_{t_{j}}^{K}=\sum_{k=0}^{M}\beta_{k}\phi_{k}
