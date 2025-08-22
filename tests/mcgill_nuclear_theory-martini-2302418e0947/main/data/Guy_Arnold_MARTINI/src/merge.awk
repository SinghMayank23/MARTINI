#!/usr/bin/awk -f
BEGIN{trlen=1;runs=ARGC-1;count=1;len=0;}
{ 
if($1!="E" && NF!=0)
    {
        trlen++;
	t[count] = $1;
        omega[count] = $2;
	pq[count] = $3;
        pg[count] = $4;
	count++;
    }
}

END{
#  printf("%i\n",count);
  for(i=1;i<=count/2;i++) 
    {
        printf("%e %e %e %e %e %e\n",t[i],omega[i],pq[i],pg[i],pq[i+(count+1)/2-1],pg[i+(count+1)/2-1]);
	if((i)%399==0 && i!=1)printf("\n");
    }
 
}

