/*************************************************************
* PROGRAMA: gi-mle-mc-v3.ox
*
* AUTORA: Gabrielle Carvalho
*
* DATA: 07/07/2022
*
* SOBRE: Simulação de Monte Carlo, para a distribuiçaõ Gaussiana
* Inversa. Para valores de T = (25,50,75,100,250,400)
*
*
*************************************************************/


#include <oxstd.oxh>
#include <oxprob.oxh>
#import <maximize>
#include <oxdraw.oxh>
#include <oxfloat.oxh>

static decl s_vy; // vetor da amostra global


//função Log-Verossimilhança
fLogLik(const vP, const adFunc, const avScore, const amHess)
{
	decl N = rows(s_vy);  // número de linhas do vetor y
	decl a = vP[0]; // mi
	decl b = vP[1]; // lambda


	 	//função log-verossimilhança
	 adFunc[0] = sumc( 0.5*log(vP[1]) - 1.5*log(s_vy)   	- (vP[1]*(s_vy - vP[0]) .^2) ./ (2*vP[0]^2*s_vy) );

		//gradiente analítico	
	if(avScore)
{

	    (avScore[0])[0] =	 sumc( (b/(a^3*s_vy)) * (a*(s_vy-a) + (s_vy - a).^2 )  );
        (avScore[0])[1] =	 sumc( (1/(2*b))  - ( ( (s_vy - a).^2 )./(2* a^2 *s_vy))  );

}	 
	if( isnan(adFunc[0]) || isdotinf(adFunc[0]) )
	return 0;
	else

	return 1; // 1 indicando sucesso
}

main()
{  	decl nobs = 400; // tamanho da amostra
	decl nrep = 10000; // número de réplicas
	decl vp, dfunc; // vetores para guardar valores dos parâmetros e da função BFGS
	decl a,b; // valores verdadeiros dos parâmetros
	decl temp;
	decl ir, i; // vetores para maximização da BFGS e para realização do laço
	decl vmlea, vmleb; // vetores para armazenar os parâmetros estimados para cada laço i
	decl hessiano, var, erro; // matriz hessiana, variância e erros-padrão
	decl cFailure = 0; // número de falhas de convergência
	decl exectime; // tempo de execução
	
	exectime = timer(); // inicia contador do tempo de execução

	ranseed("MWC_52"); // gerador de números aleatórios
	ranseed(1987); // semente
	
	/* Parâmetros verdadeiros */
	 a = 1.5;		//mi
	 b = 2;			//lambda

	/* Vetores de zeros que armazenarão valores estimador */
	vmlea = zeros(nrep, 1);
	vmleb = zeros(nrep, 1);

	/* Limitador de iterações */
	MaxControl(50, -1); // limit iterations

	/* Laço de Monte Carlo */
		 	for(i=0; i<nrep; i++)
	{
	
		s_vy  = raninvgaussian(nobs, 1, a, b); // amostra aleatória gaussiana inversa
		vp = <1.2; 1.7>; // palpite inicial

		/* Maximização via BFGS */
		//	MaxControl(-1, 1); 
		//ir = MaxBFGS(fLogLik, &vp, &dfunc, 0, FALSE);
		ir = MaxBFGS(fLogLik, &vp, &dfunc, 0, TRUE);
			if(ir == MAX_CONV || ir == MAX_WEAK_CONV)
			{
				vmlea[i] = vp[0];
				vmleb[i] = vp[1];

				/* Hessiana e erro-padrão assintótico */
				Num2Derivative(fLogLik, vp, &hessiano);
				var =  diagonal(invertsym(-hessiano))';
				erro = sqrt(diagonal(invertsym(-hessiano)))';
			}
			
			else
			{
				i--;
				cFailure++;
			}
	}

   
	print("\narmazenamento", "%6.3f ", vmlea[0]);
	print("\nAutor do Programa: Gabrielle Carvalho");
	print("\nNome do arquivo fonte: gi-mle-mc-v3.ox");
	print("\nGerador de números aleatórios: ", "MWC_52");
	print("\nData de execução: ", date());
	print("\nHorário de execução: ",time());
	print("\n");
	print("\nDistribuição: Gaussiana Inversa");
	print("\nMétodo de estimação: máxima verossimilhança");	 

	print("\n");
	
	// cálculo de médias, variâncias, desvios-padrão, viés, viés relativo (%),
    // erros quadrados médios, assimetria e curtose das estimativas de mi
	print("\nMomentos MLE para mi: ", "%10.4f", "%c", {"média  ", " desvio-padrão ",
	"assimetria ", "curtose"}, moments(vmlea)[1:4][]');
	print("\nTamanho da amostra: ", "%6d", nobs);
	print("\nNum. Rep. Monte Carlo: ", "%6d", nrep);
	print("\nNum. falhas: ", "%6d", cFailure);
	print("\nValor verdadeiro de mi: ", "%6.3f", a);
	print("\nMédia EMVs de mi: ", "%6.3f", double(meanc(vmlea)));
	print("\nVariância assintótica EMV de mi: ", "%6.3f", var[0]);
	print("\nErro Padrão Assintótico EMV de mi: ", "%6.3f", erro[0]);	
	print("\nViés de mi: ", "%6.3f", double(meanc(vmlea)-a));
	print("\nViés relativo de mi (%): ", "%6.3f", (double(meanc(vmlea))-a)/a*100);
	print("\nEQM de mi: ", "%6.3f", double(varc(vmlea)+ (meanc(vmlea)-a)^2));
	print("\nMáxima EMV de mi: ", "%6.3f", double(max(vmlea)));
	print("\nMínima EMV de mi: ", "%6.3f", double(min(vmlea)));
	print("\n");

	// cálculo de médias, variâncias, desvios-padrão, viés, viés relativo (%),
    // erros quadrados médios, assimetria e curtose das estimativas de lambda
	print("\nMomentos MLE para lambda: ", "%10.4f", "%c", {"média  ", " desvio-padrão ",
	"assimetria ", "curtose "}, moments(vmleb)[1:4][]');
	print("\nTamanho da amostra: ", "%6d", nobs);
	print("\nNum. Rep. Monte Carlo: ", "%6d", nrep);
	print("\nNum. falhas: ", "%6d", cFailure);
	print("\nValor verdadeiro de lambda: ", "%6.3f", b);
	print("\nMédia EMVs de lambda: ", "%6.3f", double(meanc(vmleb)));
	print("\nVariância assintótica EMV de lambda: ", "%6.3f", var[1]);
	print("\nErro Padrão Assintótico EMV de lambda: ", "%6.3f", erro[1]);	
	print("\nViés de lambda: ", "%6.3f", double(meanc(vmleb)-b));
	print("\nViés relativo de lambda (%): ", "%6.3f", (double(meanc(vmleb))-b)/b*100);
	print("\nEQM de lambda: ", "%6.3f", double(varc(vmleb)+ (meanc(vmleb)-b)^2));
	print("\nMáxima EMV de lambda: ", "%6.3f", double(max(vmleb)));
	print("\nMínima EMV de lambda: ", "%6.3f", double(min(vmleb)));
	print("\n");
	
	println( "\nEXECUTION TIME: ", timespan(exectime) );

	  // histograma de lambda
	DrawDensity(0, vmleb', "EMV", TRUE, TRUE);
    SaveDrawWindow("plot_l2.pdf");


		// histograma de mi
	DrawDensity(0, vmlea', "EMV", TRUE, TRUE);
    SaveDrawWindow("plot_m2.pdf");

	
	}