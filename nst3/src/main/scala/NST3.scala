package nst3

import java.lang.Math._
import Jama.Matrix

object NST3 {

  val M: Int = 30
  val alpha: Double =  0.0
  val vInf: (Double, Double) = (cos(alpha), sin(alpha))
  val delta = 0.01
  val gamma0 = 0.0
  val p = 3

  val points: Map[Int, (Double, Double)] = {
    val n1 = M / 3
    val n2 = M - n1 + 1
    val list: List[(Double, Double)] =
      (1 to n1).map(i => (1.5 + (i - 1) * 0.5 / (n1 - 1), 1.5 + (i - 1) * 0.5 / (n1 - 1))).toList :::
        (2 to n2).map(i => (2.0, 2.0 - 1.0 * (i - 1) / (n2 - 1))).toList
    list.zipWithIndex.map(e => e._2 + 1 -> e._1).toMap
  }

   val colPoints: Map[Int, (Double, Double)] = (1 until M).map(i => i -> ((points(i)._1 + points(i + 1)._1) / 2, (points(i)._2 + points(i + 1)._2) / 2)).toMap

   val normal: Map[Int, (Double, Double)] = (1 until M).map(k => k -> {
     val deltaK = sqrt((points(k + 1)._1 - points(k)._1) * (points(k + 1)._1 - points(k)._1) + (points(k + 1)._2 - points(k)._2) * (points(k + 1)._2 - points(k)._2))
     (-(points(k + 1)._2 - points(k)._2) / deltaK, (points(k + 1)._1 - points(k)._1) / deltaK)
   }).toMap

  var pointsP: Map[Int, List[(Double, Double)]] = Map(1 -> List(points(1)), 2 -> List(points(M / 3)), 3 -> List(points(M)))

  val R: (Double, Double, Double, Double) => Double = (x: Double, y: Double, xj: Double, yj: Double) => max(delta, sqrt((x - xj) * (x - xj) + (y - yj) * (y - yj)))

  val V: (Double, Double, Double, Double) => (Double, Double) = (x: Double, y: Double, xj: Double, yj: Double) => ((yj - y) / (2 * PI * R(x, y, xj, yj) * R(x, y, xj, yj)), (x - xj) / (2 * PI * R(x, y, xj, yj) * R(x, y, xj, yj)))

  var gamma: Map[Int, Double] = Map[Int, Double]()
  var smallGamma: Map[Int, Map[Int, Double]] = Map[Int, Map[Int, Double]]()

  val v: Int => (Double, Double) => (Double, Double) = (t: Int) => (x: Double, y: Double) =>
    (vInf._1 + (1 to M).map(i => gamma(i) * V(x, y, points(i)._1, points(i)._2)._1).sum + (1 to p).map(p => (1 to t).map(t => mul(smallGamma(p)(t), V(x, y, pointsP(p).apply(t)._1, pointsP(p).apply(t)._2))._1).sum).sum,
      vInf._2 + (1 to M).map(i => gamma(i) * V(x, y, points(i)._1, points(i)._2)._2).sum + (1 to p).map(p => (1 to t).map(t => mul(smallGamma(p)(t), V(x, y, pointsP(p).apply(t)._1, pointsP(p).apply(t)._2))._2).sum).sum)


  def gammaCalc: Int => Map[Int, Double] = (t: Int) => {
    val b: List[Double] = colPoints.toList.sortBy(_._1).map(ixy => -mul(vInf, normal(ixy._1)) - (for {
      p <- 1 to p
      t <- 1 to t
    } yield smallGamma(p)(t) * mul(V(ixy._2._1, ixy._2._2, pointsP(p).apply(t)._1, pointsP(p).apply(t)._2), normal(ixy._1))).sum
    ) ::: List(gamma0 - (1 to p).map(p => (1 to t).map(t => smallGamma(p)(t)).sum).sum)
    val a: Map[(Int, Int), Double] = (for {
      k <- 1 until M
      j <- 1 to M
    } yield (k, j) -> mul(V(colPoints(k)._1, colPoints(k)._2, points(j)._1, points(j)._2), normal(k))).toMap ++
      (for {
        j <- 1 to M
      } yield (M, j) -> 1.0).toMap
    val B = new Matrix(b.toArray, b.size)
    val A = {
      val array: Array[Array[Double]] = (1 to M).map(k => {
        (1 to M).map(j => a(k, j)).toArray
      }).toArray
      new Matrix(array)
    }
    A.solve(B).getRowPackedCopy.zipWithIndex.map(e => e._2 + 1 -> e._1).toMap
  }

  val dt = 0.25

  def findVP(point: (Double, Double), t: Int, index: Int): (Double, Double) = {

    val ddt = dt * dt /  Graph.r(v(t)(point._1, point._2)._1, v(t)(point._1, point._2)._2)
    val xx = point._1 + v(t)(point._1, point._2)._1*ddt
    val yy = point._2 + v(t)(point._1, point._2)._2*ddt
    (xx, yy)
  }

  val j = 25
//  val j = 75
//  val j = 100

  (0 to j).foreach(i => {
    println(i)
    gamma = gammaCalc(i)

    smallGamma = Map(
      1 -> smallGamma.getOrElse(1,Map[Int, Double]()).map(v => v._1 + 1 -> v._2).+(1 -> gamma(1)),
      2 -> smallGamma.getOrElse(2,Map[Int, Double]()).map(v => v._1 + 1 -> v._2).+(1 -> gamma(M/3)),
      3 -> smallGamma.getOrElse(3,Map[Int, Double]()).map(v => v._1 + 1 -> v._2).+(1 -> gamma(M))
    )

    pointsP = Map(
      1 -> (pointsP(1) ::: List(findVP(pointsP(1).last, i, 1))),
      2 -> (pointsP(2) ::: List(findVP(pointsP(2).last, i, 2))),
      3 -> (pointsP(3) ::: List(findVP(pointsP(3).last, i, 3)))
    )
  })

  def mul(a: (Double, Double), b: (Double, Double)) = a._1 * b._1 + a._2 * b._2
  def mul(a: Double, b: (Double, Double)) = (a * b._1, a * b._2)
  def sum(a: (Double, Double), b: (Double, Double)) = (a._1 + b._1, a._2 + b._2)


  val (x1, x2, y1, y2) = (0.05, 2.95, 0.05, 2.95)
  val h: Double = 0.05
  val n: Int = ((x2 - x1) / h).toInt

  val allPoints: List[(Double, Double)] = (for {
    x <- 0 to n
    y <- 0 to n
  } yield (x1 + x * h, y1 + y * h)).toList

  val Cp: Int => (Double, Double) => Double = (t: Int) => (x: Double, y: Double) => {

    val deltaT = 0.25

    val gammaP: Int => Double = (j: Int) => (gamma(j) - gamma(j)) / deltaT
    val smallGammaP: Int => Double = (p: Int) => smallGamma(p)(t) / deltaT
    val qp: Int => Double = (p: Int) => smallGammaP(p)
    val qj: Int => Double = (j: Int) => {
      if (j == 1) smallGammaP(1) + gammaP(1) else
        if (j == M / 3) smallGammaP(2) + gammaP(M / 2) else
        if (j == M) smallGammaP(3) + gammaP(M) else gammaP(j)
    }
    val Qm: Int => Double = (j: Int) => (1 to j).map(k => qj(k)).sum

    val xj: Int => Double = (j: Int) => (points(j + 1)._1 + points(j)._1) / 2
    val yj: Int => Double = (j: Int) => (points(j + 1)._2 + points(j)._2) / 2

    val xp: Int => Double = (p: Int) => {
      val index = if (p == 1) 1 else if (p == 2) M / 3 else M
      (pointsP(p).apply(t)._1 + points(index)._1) / 2
    }

    val yp: Int => Double = (p: Int) => {
      val index = if (p == 1) 1 else if (p == 2) M / 3 else M
      (pointsP(p).apply(t)._2 + points(index)._2) / 2
    }

    val W: (Double, Double) => (Double, Double) = (x: Double, y: Double) =>
      (vInf._1 + (1 to M).map(j => gamma(j) * V(x, y, points(j)._1, points(j)._2)._1).sum +
        (for {
          p <- 1 to p
          i <- 1 to t
        } yield smallGamma(p)(i) * V(x, y, pointsP(p).apply(i)._1, pointsP(p).apply(i)._2)._1).sum,
        vInf._2 + (1 to M).map(j => gamma(j) * V(x, y, points(j)._1, points(j)._2)._2).sum +
          (for {
            p <- 1 to p
            i <- 1 to t
          } yield smallGamma(p)(i) * V(x, y, pointsP(p).apply(i)._1, pointsP(p).apply(i)._2)._2).sum)

    def dphi(): Double = {
      val s1 = (1 to M - 1).map(j => {
        Qm(j) * ((yj(j) - y)*(points(j + 1)._1 - points(j)._1) + (x - xj(j))*(points(j + 1)._2 - points(j)._2)) / (2 * PI * ((x - xj(j)) * (x - xj(j)) + (y - yj(j)) * (y - yj(j))))
      }).sum

      val s2 = (1 to p).map(p => {
        val index = if (p == 1) 1 else if (p == 2) M / 3 else M
        qp(p)* ((yp(p) - y) * (points(index)._1 - pointsP(p).apply(t)._1) + (x - xp(p)) * (points(index)._2 - pointsP(p).apply(t)._2)) /
          (2 * PI * ((x - xp(p)) * (x - xp(p)) + (y - yp(p)) * (y - yp(p))))
      }).sum

      val s3 = (for {
        p <- 1 to p
        i <- 1 to t - 1
      } yield {
        smallGamma(p)(i)* mul(V(x, y, pointsP(p).apply(i)._1, pointsP(p).apply(i)._2), W(pointsP(p).apply(i)._1, pointsP(p).apply(i)._2))
      }).sum

      s1 + s2 - s3
    }

    val cp = 1 - mul(v(t)(x, y), v(t)(x, y)) / mul(vInf, vInf) - (2 / mul(vInf, vInf)) * dphi()

    cp
  }

  def main(args: Array[String]): Unit = {
    Graph.paintPlot(Cp(j), allPoints, points.values.toList, pointsP.values.flatten.toList, "CP - V", Some(v(j)))
  }
}
