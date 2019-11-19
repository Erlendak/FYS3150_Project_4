# FYS3150_Project_4

I dette projektet har vi tatt for oss en elementær model for å beskrive termodynamikk, i et todimensjonalt binært system kalt isingmodelen. Gjennom Monte Carlo integrasjon og ved hjelp av metropolis algoritmen har vi simulert utviklingen av system med forskjelige størrelser for å kunne approksimere systemet når størrelsen går mot uendelig, i forskjellige temperaturer.

Dette har vært et interessant og spenndene projekt, men også et krevende projekt. Ettersom at beregningene i prosjektet er ganske tunge, og vi har hatt et behov for å parallellisere koden. Siden vi begge hovedsakelig bruker Windows valgte vi å implementere parallelliseringen gjennom Open MP, noe som ikke ble velykket, og vi endte opp med å simulere i serie i stenden.

Resultatene våre tyder på at vi trenger omtrent $10^4$ eller mer Monte Carlo iterasjoner for å oppnå gode resultater, for modellene vi foretar oss i projectet. Ettersom modellene da oppnår et equilibrium stadie. I tillegg viser resultatene våre at høyere temperaturer aksepterer flere forskjellige konfigurasjoner, og dermed hopper lettere mellom forskjellige energier, enn kalde temperaturer som har en tendens til å holde én spesefikk energi, der alle spin peker i én rettning. Med andre ord har vi funnet ut at modellen følger en Boltzmanns distribusjon. 
