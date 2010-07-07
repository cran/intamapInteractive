library(intamapInteractive)

data(walker)
coordinates(walker)=~X+Y
object=createIntamapObject(observations=walker)
object=anisotropyChoice(object)

print(summary(object$clusters$index))
print(object$anisPar)



object=doSegmentation(object)

print(summary(object$clusters$index))