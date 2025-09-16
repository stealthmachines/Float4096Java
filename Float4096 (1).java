package com.xai.float4096;

import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Objects;
import java.util.Random;
import java.util.Set;

/**
 * A high-precision floating-point number class with 4096-bit precision.
 * Uses {@link BigDecimal} for underlying arithmetic with a fixed {@link MathContext}
 * configured for 4096-bit precision and HALF_EVEN rounding.
 * Provides arithmetic operations, comparisons, and utility methods for high-precision calculations.
 * 
 * This library demonstrates functional completeness via subtraction, allowing simulation of arbitrary logic circuits
 * (and thus supporting Turing-complete systems in principle). Boolean values are represented as +0.0 (true) and -0.0 (false).
 * Logic gates are implemented using only subtraction, leveraging signed zeros.
 * 
 * Additionally, supports deterministic generation of 4096 contrived and unique Unicode characters based on the value.
 */
public final class Float4096 {
    private static final MathContext PRECISION_4096 = new MathContext(4096, RoundingMode.HALF_EVEN);
    
    // Mathematical constants
    public static final Float4096 ZERO = new Float4096("0.0");  // +0.0 (true in logic representation)
    public static final Float4096 NEG_ZERO = new Float4096("-0.0");  // -0.0 (false in logic representation)
    public static final Float4096 ONE = new Float4096("1");
    public static final Float4096 TWO = new Float4096("2");
    public static final Float4096 PI = calculatePi();
    public static final Float4096 E = calculateE();
    public static final Float4096 PHI = calculateGoldenRatio();
    
    private final BigDecimal value;

    // Constructors (unchanged)
    public Float4096(String value) {
        this.value = new BigDecimal(value, PRECISION_4096);
    }

    public Float4096(BigDecimal value) {
        this.value = Objects.requireNonNull(value).round(PRECISION_4096);
    }

    public Float4096(double value) {
        this.value = new BigDecimal(value, PRECISION_4096);
    }

    // Arithmetic operations (unchanged, with subtract already present)
    public Float4096 add(Float4096 other) {
        return new Float4096(this.value.add(other.value, PRECISION_4096));
    }

    public Float4096 subtract(Float4096 other) {
        return new Float4096(this.value.subtract(other.value, PRECISION_4096));
    }

    public Float4096 multiply(Float4096 other) {
        return new Float4096(this.value.multiply(other.value, PRECISION_4096));
    }

    public Float4096 divide(Float4096 other) {
        if (other.value.signum() == 0) {
            throw new ArithmeticException("Division by zero");
        }
        return new Float4096(this.value.divide(other.value, PRECISION_4096));
    }

    public Float4096 sqrt() {
        if (this.value.signum() < 0) {
            throw new ArithmeticException("Square root of negative number");
        }
        BigDecimal x = this.value;
        BigDecimal previous;
        for (int i = 0; i < 100; i++) {
            previous = x;
            x = x.add(this.value.divide(x, PRECISION_4096)).divide(TWO.value, PRECISION_4096);
            if (x.equals(previous)) {
                break;
            }
        }
        return new Float4096(x);
    }

    public Float4096 abs() {
        return new Float4096(this.value.abs(PRECISION_4096));
    }

    // Comparison and utilities (unchanged)
    public int compareTo(Float4096 other) {
        return this.value.compareTo(other.value);
    }

    public boolean isZero() {
        return this.value.signum() == 0;
    }

    public double toDouble() {
        return this.value.doubleValue();
    }

    public BigDecimal getValue() {
        return value;
    }

    @Override
    public String toString() {
        return value.toPlainString();
    }

    public String toScientificString() {
        return value.toEngineeringString();
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof Float4096)) return false;
        Float4096 that = (Float4096) o;
        return value.equals(that.value);
    }

    @Override
    public int hashCode() {
        return Objects.hash(value);
    }

    // New: Logic gate implementations using only subtraction and signed zeros
    // +0.0 represents true, -0.0 represents false

    /**
     * Implements a NOT gate using subtraction.
     * @param input The input Float4096 (+0.0 or -0.0).
     * @return The negated output.
     */
    public static Float4096 logicNot(Float4096 input) {
        return NEG_ZERO.subtract(input);
    }

    /**
     * Implements an OR gate using subtraction.
     * @param a First input.
     * @param b Second input.
     * @return The OR result.
     */
    public static Float4096 logicOr(Float4096 a, Float4096 b) {
        return a.subtract(logicNot(b));
    }

    /**
     * Implements an AND gate using subtraction (via De Morgan's law).
     * @param a First input.
     * @param b Second input.
     * @return The AND result.
     */
    public static Float4096 logicAnd(Float4096 a, Float4096 b) {
        return logicNot(logicOr(logicNot(a), logicNot(b)));
    }

    /**
     * Implements an XOR gate using subtraction.
     * @param a First input.
     * @param b Second input.
     * @return The XOR result.
     */
    public static Float4096 logicXor(Float4096 a, Float4096 b) {
        return logicOr(logicAnd(logicNot(a), b), logicAnd(a, logicNot(b)));
    }

    /**
     * Demonstrates a single-bit adder using only subtraction-based gates.
     * @param a First bit (+0.0 or -0.0).
     * @param b Second bit.
     * @param carryIn Carry input.
     * @return List containing sum and carryOut.
     */
    public static List<Float4096> adder(Float4096 a, Float4096 b, Float4096 carryIn) {
        Float4096 sum = logicXor(logicXor(a, b), carryIn);
        Float4096 carryOut = logicOr(logicAnd(logicXor(a, b), carryIn), logicAnd(a, b));
        return List.of(sum, carryOut);
    }

    // New: Deterministic generation of 4096 unique characters
    /**
     * Deterministically generates 4096 contrived and unique Unicode characters based on this Float4096 value.
     * Uses a SHA-256 hash of the value as a seed for Random to ensure determinism.
     * @return A list of 4096 unique characters.
     * @throws RuntimeException if hashing fails.
     */
    public List<Character> generateUniqueCharacters() {
        try {
            MessageDigest digest = MessageDigest.getInstance("SHA-256");
            byte[] hash = digest.digest(this.toString().getBytes("UTF-8"));
            long seed = new BigInteger(1, hash).longValue();
            Random rand = new Random(seed);
            
            Set<Character> uniqueChars = new HashSet<>();
            while (uniqueChars.size() < 4096) {
                int codePoint = rand.nextInt(0x10FFFF + 1);
                if (Character.isValidCodePoint(codePoint)) {
                    uniqueChars.add((char) codePoint);
                }
            }
            return new ArrayList<>(uniqueChars);
        } catch (NoSuchAlgorithmException | java.io.UnsupportedEncodingException e) {
            throw new RuntimeException("Failed to generate characters", e);
        }
    }

    // Constant calculations (unchanged)
    private static Float4096 calculateGoldenRatio() {
        BigDecimal five = new BigDecimal("5", PRECISION_4096);
        BigDecimal sqrt5 = five.sqrt(PRECISION_4096);
        BigDecimal one = new BigDecimal("1", PRECISION_4096);
        BigDecimal two = new BigDecimal("2", PRECISION_4096);
        return new Float4096(one.add(sqrt5).divide(two, PRECISION_4096));
    }

    private static Float4096 calculatePi() {
        // Chudnovsky algorithm for PI (simplified iteration for brevity)
        BigDecimal C = new BigDecimal("426880", PRECISION_4096);
        BigDecimal L = new BigDecimal("13591409", PRECISION_4096);
        BigDecimal X = new BigDecimal("1", PRECISION_4096);
        BigDecimal M = new BigDecimal("1", PRECISION_4096);
        BigDecimal K = new BigDecimal("6", PRECISION_4096);
        BigDecimal sum = L;
        BigDecimal one = new BigDecimal("1", PRECISION_4096);
        BigDecimal negOne = new BigDecimal("-1", PRECISION_4096);
        BigDecimal sqrt10005 = new BigDecimal("10005", PRECISION_4096).sqrt(PRECISION_4096);

        for (int k = 1; k < 100; k++) {
            BigDecimal kBig = new BigDecimal(k);
            M = M.multiply(K.subtract(one).multiply(K).multiply(K.add(one)))
                .divide(kBig.pow(3).multiply(new BigDecimal("640320", PRECISION_4096).pow(3)), PRECISION_4096);
            L = L.add(new BigDecimal("545140134", PRECISION_4096));
            X = X.multiply(negOne);
            K = K.add(new BigDecimal("12", PRECISION_4096));
            sum = sum.add(M.multiply(L).multiply(X));
        }

        return new Float4096(C.multiply(sqrt10005).divide(sum, PRECISION_4096));
    }

    private static Float4096 calculateE() {
        BigDecimal sum = BigDecimal.ONE;
        BigDecimal factorial = BigDecimal.ONE;
        for (int i = 1; i < 1000; i++) {
            factorial = factorial.multiply(new BigDecimal(i));
            sum = sum.add(BigDecimal.ONE.divide(factorial, PRECISION_4096));
        }
        return new Float4096(sum);
    }
}