import chisel3._
import chisel3.util._

class MulIO(val len: Int) extends Bundle {
  val in = Flipped(ValidIO(Vec(2, Output(SInt(len.W)))))
  val flush  = Input(Bool())
  val out = ValidIO(Output(SInt((len * 2).W)))
}

class FullAdder(n: Int = 32) extends Module {
  val io = IO(new Bundle {
    val a = Input(SInt(n.W))
    val b = Input(SInt(n.W))
    val cin = Input(SInt(n.W))
    val sum = Output(SInt(n.W))
    val cout = Output(SInt(n.W))
  })

  // Generate the sum
  val a_xor_b = io.a ^ io.b
  io.sum := a_xor_b ^ io.cin
  // Generate the carry
  val a_and_b = io.a & io.b
  val b_and_cin = io.b & io.cin
  val a_and_cin = io.a & io.cin
  io.cout := a_and_b | b_and_cin | a_and_cin
}


class Booth4WTMultiplier(slen: Int = 33) extends Module {
  val io = IO(new MulIO(len = slen))
  val len = slen - 1
  val fa  = Seq.fill(len/4 - 1)(Module(new FullAdder(2*slen)))

  val (x, y) = (io.in.bits(0), io.in.bits(1))

  val px = WireInit(VecInit(Seq.fill(len/4)(0.S((2*slen).W))))
  val prod = WireInit(VecInit(Seq.fill(len/4 + 1)(0.S((2*slen).W))))
  val y_ext = Cat(y.asSInt(),0.S(1.W))
  val y_comp = WireInit(VecInit(Seq.fill(len/4)(0.U((len/4).W))))

  val m8condi = WireInit(VecInit(Seq.fill(len/4)(false.B)))

  val mx0 = WireInit(VecInit(Seq.fill(len/4)(0.S((slen+3).W))))
  val mx1 = WireInit(VecInit(Seq.fill(len/4)(0.S((slen+3).W))))

  val x1condi = WireInit(VecInit(Seq.fill(len/4)(false.B)))
  val x2condi = WireInit(VecInit(Seq.fill(len/4)(false.B)))
  val x3condi = WireInit(VecInit(Seq.fill(len/4)(false.B)))
  val x2 = (x << 1.U).asSInt()
  val x3 = x + x2
  val xminus = WireInit(VecInit(Seq.fill(len/4)(0.S((slen+2).W))))
  for(j <- 0 until (len/4)) {
    val i = 4 * j + 1

    y_comp(j) := Mux(y_ext(i+3), (~y_ext(i+2, i-1)).asUInt(), y_ext(i+2, i-1).asUInt())

    m8condi(j) := (y_comp(j)(0) || y_comp(j)(1) || y_comp(j)(2)) && y_comp(j)(3)
    mx0(j) := x << (2.U(2.W) + m8condi(j).asUInt())

    x1condi(j) := y_comp(j)(0) ^ y_comp(j)(1)
    x2condi(j) := (y_comp(j)(0) || y_comp(j)(1)) ^ y_comp(j)(2)

    when(x1condi(j) && x2condi(j)) {
      xminus(j) := x3
    }.elsewhen(!x1condi(j) && x2condi(j)) {
      xminus(j) := x2
    }.elsewhen(x1condi(j) && !x2condi(j)){
      xminus(j) := x
    }.otherwise {
      xminus(j) := 0.S((slen+2).W)
    }

    mx1(j) := Mux(y_comp(j).orR, mx0(j) - xminus(j), 0.S((slen+3).W))

    printf("mx0 %d\n", mx0(j))
    printf("mx1 %d\n", mx1(j))

    px(j) := Mux(y_ext(i+3), -mx1(j), mx1(j))
    printf("px %d\n", px(j))
  }
  val mux1 = Mux(y_ext(32)^y_ext(33),x,0.S((slen).W))

  prod(0) := (px(0)).asSInt()
  prod(1) := (px(1)<<4).asSInt()
  prod(2) := (px(2)<<8).asSInt()
  prod(3) := (px(3)<<12).asSInt()
  prod(4) := (px(4)<<16).asSInt()
  prod(5) := (px(5)<<20).asSInt()
  prod(6) := (px(6)<<24).asSInt()
  prod(7) := (px(7)<<28).asSInt()
  prod(8) := (Mux(y_ext(32), mux1, -mux1)<<32).asSInt()

  //----------------------------WT Structure---------------------------
  //branch0
  fa(0).io.a   := prod(0)
  fa(0).io.b   := prod(1)
  fa(0).io.cin := prod(2)
  fa(1).io.a   := prod(3)
  fa(1).io.b   := prod(4)
  fa(1).io.cin := prod(5)
  fa(2).io.a   := prod(6)
  fa(2).io.b   := prod(7)
  fa(2).io.cin := prod(8)

  //branch1
  fa(3).io.a   := fa(0).io.sum
  fa(3).io.b   := fa(1).io.sum
  fa(3).io.cin := fa(2).io.sum
  fa(4).io.a   := fa(0).io.cout<<1
  fa(4).io.b   := fa(1).io.cout<<1
  fa(4).io.cin := fa(2).io.cout<<1

  //branch2
  fa(5).io.a   := fa(3).io.sum
  fa(5).io.b   := fa(4).io.sum
  fa(5).io.cin := fa(3).io.cout<<1

  //branch3
  fa(6).io.a   := fa(5).io.sum
  fa(6).io.b   := fa(4).io.cout<<1
  fa(6).io.cin := fa(5).io.cout<<1

  printf("fa(6).io.sum %d\n",fa(6).io.sum)
  printf("fa(6).io.cout %d\n",fa(6).io.cout)
  printf("fa(6).io.cout<<1 %d\n",fa(6).io.cout<<1)
  val sumpx = fa(6).io.sum + (fa(6).io.cout<<1)

  //   val sumpx = prod(0) + prod(1) + prod(2) + prod(3) +
  //     prod(4) + prod(5) + prod(6) + prod(7) + prod(8)

  printf("sumpx %d", sumpx)
  io.out.bits :=  sumpx
  io.out.valid := io.in.valid
}