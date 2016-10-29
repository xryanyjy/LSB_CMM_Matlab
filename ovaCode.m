function ova_value = ovaCode(L,K)

code = zeros(L,K);

for i = 0:(K-1)
    a(i+1) = i;
end
cls = mod(a,L)+1;

for k = 1:K
    code(cls(k),k) = 1;
end

code = code .* 10 + 1;
ova_value.ealpha = code;
ova_value.cls = cls;


