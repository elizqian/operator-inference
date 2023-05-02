function D2 = dupmat(n)
  m   = n * (n + 1) / 2;
  nsq = n^2;
  r   = 1;
  a   = 1;
  v   = zeros(1, nsq);
  cn  = cumsum(n:-1:2);   % [EDITED, 2021-08-04], 10% faster
  for i = 1:n
     % v(r:r + i - 2) = i - n + cumsum(n - (0:i-2));
     v(r:r + i - 2) = i - n + cn(1:i - 1);   % [EDITED, 2021-08-04]
     r = r + i - 1;
     
     v(r:r + n - i) = a:a + n - i;
     r = r + n - i + 1;
     a = a + n - i + 1;
  end
  
  D2 = sparse(1:nsq, v, 1, nsq, m);
end