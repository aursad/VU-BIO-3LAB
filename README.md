VU-BIO-3LAB
===========
<b>Užduotis</b><br/>
===========
Sukurti pobes parinkimo diagnostinei sistemai programą, skirtą sukurti probe tinkama tik tam tikrų virusų DNR detekcijai, kuomet padauginimui naudojma PGR (Polimerazinę grandininę reakciją)
Sukurtą programą pritaikyti projektuojant diagnostinę sistemą žmogaus papilomos viruso (HPV) 16,18,31,33,35,51,52 tipų (sukelia gimdos kaklelio vėžį) infekcijai nustatyti. Diagnostinė sistema turi būti specifiška ir turi nebūti signalo nepavojingų HPV infekcijų atvejais (tipai 6, 11, 40, 42, 43, 44, 57, 81).
Dokumentuota programa turėtų parinkinėti DNR probe L1 geno regione. Atitinkama šio geno seka iš HPV 16 tipo genomo pateikta šio dokumento pabaigoje. Šią seką galite naudoti per blast užklausas rinkdami sekas analizei.
Kuriant algoritmą tikslas būtų sukurti jį tokį, kad būtų parenkamas kuo mažesnis probiu kiekis tinkantis visų žinomų didelės rizikos HPV tipų variantų diagnostikai.
<br/><br/>

I dalies uzduotis:
Modifikuokite jau sukurta  programa (Ankstesnis laboratorinis) taip, kad:
<ol>
<li>
Parsiustu visu tipu pavojingu ir nepavojingu zmogaus papilomos viruso prieinamas sekas fasta formatu į atskirus failus (vieno tipo sekos vienas failas). Paieskai naudokite dokumento pabaigoje pateikta L1 geno fragmenta.
</li><li>
Naudojant cd-hit programa ( http://weizhong-lab.ucsd.edu/cd-hit/) panaikintu identiskas sekas.
</li><li>
Sulietus atskirus failus į viena fasta faila (virsuje pavojingi tipai, zemiau nepavojingi) ir gautumet bendra sulyginima su mafft (http://mafft.cbrc.jp/alignment/software/)
</li></ol>
