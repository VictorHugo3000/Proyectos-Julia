"""
    fibonacci(n)

Calcula los primeros `n` términos de la sucesión de *Fibonacci* y los entrega en un vector de tamaño `n`.

## Ejemplo
```
fibonacci(13)
```
"""
function fibonacci(n) # n = número de términos de la sucesión
    if n == 1
        a = [0]
    else
        a = zeros(Int, n)
        a[2] = 1
        if n ≥ 3
            for i ∈ 3:n
                a[i] = a[i-1] + a[i-2]
            end
        end
    end
    return a
end

"""
    fibonacciMax(m)

Calcula los términos de la sucesión de *Fibonacci* hasta un valor no mayor que `m`.

## Ejemplo
fibonacciMax(145)
"""
function fibonacciMax(m) # m = valor máximo de la sucesión
    if m == 0
        return [0]
    else
        a = [0, 1]
        while a[end] ≤ m
            n = a[end] + a[end - 1]
            push!(a, n)
        end
        return a[1:(end - 1)]
    end
end

"""
    bisección(f, a, b; δ = abs((a + b)/2) / 1_000_000, m = 10_000)

Aplica el *algoritmo de bisección* para encontrar la raíz de una función real y continua `f` de variable real en un intervalo cerrado `[a,b]` bajo el supuesto de que la raíz es única en dicho itervalo y que f(a)f(b) < 0. Opcionalmente se puede especificar un nivel de tolerancia `δ` y/o un número máximo `m` de iteraciones, si no se especifican entonces se utilizan los valores por defecto que se muestran. El resultado es una **tupla** etiquetada con los valores de la raíz encontrada, el valor de la función en dicha raíz, el número de iteraciones realizadas, el número máximo de iteraciones, y la tolerancia utilizada. En caso de que el algoritmo entregue un resultado, no porque se cumplió el *criterio de paro* sino porque se alcanzó el número máximo de iteraciones, se mostrará un mensaje de advertencia al respecto.

## Ejemplo 1
```
f(x) = (x - 3) * (x - 1) * (x + 1)
bisección(f, -1.9, 0.0)
```

## Ejemplo 2
```
f(x) = (x - 3) * (x - 1) * (x + 1)
bisección(f, -1.9, 0.0, m = 9)
```
"""
function bisección(f, a, b; δ = abs((a + b)/2) / 1_000_000, m = 10_000)
    iter = 1
    z = (a + b) / 2
    while iter ≤ m
        if f(z) == 0.0 || (b - a)/2 ≤ δ
            break # criterio de paro
        elseif f(a) * f(z) > 0
            a = z
        else
            b = z
        end
        z = (a + b) / 2
        iter += 1
    end
    iter -= 1
    if iter == m
        println("Se alcanzó el número máximo de iteraciones = $m")
    end
    return (raíz = z, dif = f(z), numiter = iter, maxiter = m, tol = δ)
end

function MisFunciones()
    println("bisección", "\t\t", "Algoritmo de bisección para encontrar raíces")
    println("fibonacci", "\t\t", "Generar cierto número de términos de Fibonacci")
    println("fiboncciMax", "\t\t", "Generar términos de Fibonacci que no excedan un valor")
    println("\n", "Más detalle: Ejecuta ?nombre donde `nombre` es cualquiera de las funciones")
    return nothing
end

println("funciones: bisección   fibonacci   fibonacciMax")
println("Ejecuta MisFunciones() para resumen")