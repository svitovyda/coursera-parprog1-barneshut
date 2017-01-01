import barneshut.conctrees._
import common._

package object barneshut {

  val SECTOR_PRECISION = 8
  def minimumSize = 0.00001f
  def delta: Float = 0.01f
  def theta = 0.5f
  def eliminationThreshold = 0.5f
  def force(m1: Float, m2: Float, dist: Float): Float = gee * m1 * m2 / (dist * dist)
  def gee: Float = 100.0f
  def distance(x0: Float, y0: Float, x1: Float, y1: Float): Float = {
    math.sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0)).toFloat
  }

  sealed abstract class Quad {
    def massX: Float

    def massY: Float

    def mass: Float

    def centerX: Float

    def centerY: Float

    def size: Float

    def total: Int

    def insert(b: Body): Quad
  }

  class Boundaries {
    var minX = Float.MaxValue

    var minY = Float.MaxValue

    var maxX = Float.MinValue

    var maxY = Float.MinValue
    def size: Float = math.max(width, height)
    def centerX: Float = minX + width / 2
    def width: Float = maxX - minX
    def centerY: Float = minY + height / 2
    def height: Float = maxY - minY
    override def toString = s"Boundaries($minX, $minY, $maxX, $maxY)"
  }

  case class Empty(centerX: Float, centerY: Float, size: Float) extends Quad {
    def massX: Float = centerX
    def massY: Float = centerY
    def mass: Float = 0
    def total: Int = 0
    def insert(b: Body): Quad = Leaf(centerX, centerY, size, Seq(b))
  }

  case class Fork(
    nw: Quad, ne: Quad, sw: Quad, se: Quad
  ) extends Quad {
    val centerX: Float = (nw.centerX + ne.centerX) / 2
    val centerY: Float = (nw.centerY + sw.centerY) / 2
    val size: Float = nw.size + ne.size
    val mass: Float = List(nw, ne, sw, se).map(_.mass).sum
    val massX: Float = if (mass == 0) centerX else List(nw, ne, sw, se).map(q => q.mass * q.massX).sum / mass
    val massY: Float = if (mass == 0) centerY else List(nw, ne, sw, se).map(q => q.mass * q.massY).sum / mass
    val total: Int = List(nw, ne, sw, se).map(_.total).sum

    def insert(b: Body): Fork = {
      val dx = b.x > centerX
      val dy = b.y > centerY
      (dx, dy) match {
        case (false, false) => Fork(nw.insert(b), ne, sw, se)
        case (true, false)  => Fork(nw, ne.insert(b), sw, se)
        case (false, true)  => Fork(nw, ne, sw.insert(b), se)
        case (true, true)   => Fork(nw, ne, sw, se.insert(b))
      }
    }
  }

  case class Leaf(centerX: Float, centerY: Float, size: Float, bodies: Seq[Body]) extends Quad {
    val mass: Float = bodies.map(_.mass).sum
    val massX: Float = bodies.map(b => b.mass * b.x).sum / mass
    val massY: Float = bodies.map(b => b.mass * b.y).sum / mass
    val total: Int = bodies.length

    def insert(b: Body): Quad = {
      val next = bodies :+ b
      if (size > minimumSize) {
        val quadSize = size / 4
        val halfSize = size / 2
        val nw = Empty(centerX - quadSize, centerY - quadSize, halfSize)
        val ne = Empty(centerX + quadSize, centerY - quadSize, halfSize)
        val sw = Empty(centerX - quadSize, centerY + quadSize, halfSize)
        val se = Empty(centerX + quadSize, centerY + quadSize, halfSize)
        next.foldLeft(Fork(nw, ne, sw, se))((agg, body) => agg.insert(body))
      } else copy(bodies = next)
    }
  }

  class Body(val mass: Float, val x: Float, val y: Float, val xspeed: Float, val yspeed: Float) {

    def updated(quad: Quad): Body = {
      var netforcex = 0.0f
      var netforcey = 0.0f

      def addForce(thatMass: Float, thatMassX: Float, thatMassY: Float): Unit = {
        val dist = distance(thatMassX, thatMassY, x, y)
        /* If the distance is smaller than 1f, we enter the realm of close
         * body interactions. Since we do not model them in this simplistic
         * implementation, bodies at extreme proximities get a huge acceleration,
         * and are catapulted from each other's gravitational pull at extreme
         * velocities (something like this:
         * http://en.wikipedia.org/wiki/Interplanetary_spaceflight#Gravitational_slingshot).
         * To decrease the effect of this gravitational slingshot, as a very
         * simple approximation, we ignore gravity at extreme proximities.
         */
        if (dist > 1f) {
          val dforce = force(mass, thatMass, dist)
          val xn = (thatMassX - x) / dist
          val yn = (thatMassY - y) / dist
          val dforcex = dforce * xn
          val dforcey = dforce * yn
          netforcex += dforcex
          netforcey += dforcey
        }
      }

      def traverse(quad: Quad): Unit = (quad: Quad) match {
        case Empty(_, _, _) => // no force
        case Leaf(_, _, _, bodies) => // add force contribution of each body by calling addForce
          bodies.foreach(b => addForce(b.mass, b.x, b.y))
        case Fork(nw, ne, sw, se) => // see if node is far enough from the body, or recursion is needed
          val dist = distance(quad.massX, quad.massY, x, y)
          if (quad.size / dist < theta) addForce(quad.mass, quad.massX, quad.massY)
          else Seq(nw, ne, sw, se).foreach { q => traverse(q) }
      }

      traverse(quad)

      val nx = x + xspeed * delta
      val ny = y + yspeed * delta
      val nxspeed = xspeed + netforcex / mass * delta
      val nyspeed = yspeed + netforcey / mass * delta

      new Body(mass, nx, ny, nxspeed, nyspeed)
    }

  }

  class SectorMatrix(val boundaries: Boundaries, val sectorPrecision: Int) {
    val sectorSize: Float = boundaries.size / sectorPrecision
    val matrix = new Array[ConcBuffer[Body]](sectorPrecision * sectorPrecision)
    for (i <- matrix.indices) matrix(i) = new ConcBuffer

    def +=(b: Body): SectorMatrix = {
      val x = math.min(math.max(boundaries.minX, b.x), boundaries.maxX)
      val y = math.min(math.max(boundaries.minY, b.y), boundaries.maxY)
      apply(((x - boundaries.minX) / sectorSize).toInt, ((y - boundaries.minY) / sectorSize).toInt) += b
      this
    }

    def apply(x: Int, y: Int): ConcBuffer[Body] = matrix(y * sectorPrecision + x)

    def combine(that: SectorMatrix): SectorMatrix = {
      for (i <- matrix.indices) matrix(i) = matrix(i).combine(that.matrix(i))
      this
    }

    def toQuad(parallelism: Int): Quad = {
      def BALANCING_FACTOR = 4
      def quad(x: Int, y: Int, span: Int, achievedParallelism: Int): Quad = {
        if (span == 1) {
          val sectorSize = boundaries.size / sectorPrecision
          val centerX = boundaries.minX + x * sectorSize + sectorSize / 2
          val centerY = boundaries.minY + y * sectorSize + sectorSize / 2
          var emptyQuad: Quad = Empty(centerX, centerY, sectorSize)
          val sectorBodies = this (x, y)
          sectorBodies.foldLeft(emptyQuad)(_ insert _)
        } else {
          val nspan = span / 2
          val nAchievedParallelism = achievedParallelism * 4
          val (nw, ne, sw, se) =
            if (parallelism > 1 && achievedParallelism < parallelism * BALANCING_FACTOR) parallel(
              quad(x, y, nspan, nAchievedParallelism),
              quad(x + nspan, y, nspan, nAchievedParallelism),
              quad(x, y + nspan, nspan, nAchievedParallelism),
              quad(x + nspan, y + nspan, nspan, nAchievedParallelism)
            ) else (
              quad(x, y, nspan, nAchievedParallelism),
              quad(x + nspan, y, nspan, nAchievedParallelism),
              quad(x, y + nspan, nspan, nAchievedParallelism),
              quad(x + nspan, y + nspan, nspan, nAchievedParallelism)
            )
          Fork(nw, ne, sw, se)
        }
      }

      quad(0, 0, sectorPrecision, 1)
    }

    override def toString = s"SectorMatrix(#bodies: ${matrix.map(_.size).sum})"
  }

  class TimeStatistics {
    private val timeMap = collection.mutable.Map[String, (Double, Int)]()

    def clear(): Unit = timeMap.clear()

    def timed[T](title: String)(body: => T): T = {
      var res: T = null.asInstanceOf[T]
      val totalTime = /*measure*/ {
        val startTime = System.currentTimeMillis()
        res = body
        System.currentTimeMillis() - startTime
      }

      timeMap.get(title) match {
        case Some((total, num)) => timeMap(title) = (total + totalTime, num + 1)
        case None               => timeMap(title) = (0.0, 0)
      }

      println(s"$title: $totalTime ms; avg: ${timeMap(title)._1 / timeMap(title)._2}")
      res
    }

    override def toString: String = {
      timeMap map {
        case (k, (total, num)) => k + ": " + (total / num * 100).toInt / 100.0 + " ms"
      } mkString "\n"
    }
  }
}
