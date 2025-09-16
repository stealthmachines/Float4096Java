package com.xai.float4096;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;
import java.math.RoundingMode;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
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
 * 
 * Now includes quaternary (base-4) DNA logic: 
 * - Quaternary logic gates (MIN, MAX, Inverter, Successor) for multi-valued logic.
 * - DNA sequences as quaternary representations (A=0, C=1, G=2, T=3).
 * - Safe math operations on DNA sequences (element-wise gates, addition in base-4).
 * - Base-4096 encoding/decoding using DNA sequences, where 6 DNA bases encode one base-4096 symbol (since 4^6 = 4096).
 * - Uses a fixed alphabet of 4096 unique characters for base-4096 symbols.
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

    // Arithmetic operations (unchanged)
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

    // Binary logic gates using subtraction (unchanged)
    public static Float4096 logicNot(Float4096 input) {
        return NEG_ZERO.subtract(input);
    }

    public static Float4096 logicOr(Float4096 a, Float4096 b) {
        return a.subtract(logicNot(b));
    }

    public static Float4096 logicAnd(Float4096 a, Float4096 b) {
        return logicNot(logicOr(logicNot(a), logicNot(b)));
    }

    public static Float4096 logicXor(Float4096 a, Float4096 b) {
        return logicOr(logicAnd(logicNot(a), b), logicAnd(a, logicNot(b)));
    }

    public static List<Float4096> adder(Float4096 a, Float4096 b, Float4096 carryIn) {
        Float4096 sum = logicXor(logicXor(a, b), carryIn);
        Float4096 carryOut = logicOr(logicAnd(logicXor(a, b), carryIn), logicAnd(a, b));
        return List.of(sum, carryOut);
    }

    // Deterministic generation of 4096 unique characters (now instance method, but with static fixed alphabet for base-4096)
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

    // New: Quaternary DNA Logic Section

    // DNA bases to quaternary digit mapping (A=0, C=1, G=2, T=3)
    private static final Map<Character, Integer> DNA_TO_QUAT = new HashMap<>();
    private static final Map<Integer, Character> QUAT_TO_DNA = new HashMap<>();
    static {
        DNA_TO_QUAT.put('A', 0);
        DNA_TO_QUAT.put('C', 1);
        DNA_TO_QUAT.put('G', 2);
        DNA_TO_QUAT.put('T', 3);
        QUAT_TO_DNA.put(0, 'A');
        QUAT_TO_DNA.put(1, 'C');
        QUAT_TO_DNA.put(2, 'G');
        QUAT_TO_DNA.put(3, 'T');
    }

    // Fixed base-4096 alphabet (using PHI as seed for determinism)
    private static final List<Character> BASE4096_ALPHABET = PHI.generateUniqueCharacters();

    // Validate DNA sequence
    private static void validateDNA(String dna) {
        if (dna == null || dna.isEmpty()) {
            throw new IllegalArgumentException("Invalid DNA sequence: null or empty");
        }
        for (char c : dna.toCharArray()) {
            if (!DNA_TO_QUAT.containsKey(Character.toUpperCase(c))) {
                throw new IllegalArgumentException("Invalid DNA base: " + c);
            }
        }
    }

    // Quaternary logic gates (single values 0-3)
    /**
     * Quaternary MIN gate (analogous to AND in binary).
     * @param a First value (0-3).
     * @param b Second value (0-3).
     * @return min(a, b).
     */
    public static int quatMin(int a, int b) {
        if (a < 0 || a > 3 || b < 0 || b > 3) throw new IllegalArgumentException("Quaternary values must be 0-3");
        return Math.min(a, b);
    }

    /**
     * Quaternary MAX gate (analogous to OR in binary).
     * @param a First value (0-3).
     * @param b Second value (0-3).
     * @return max(a, b).
     */
    public static int quatMax(int a, int b) {
        if (a < 0 || a > 3 || b < 0 || b > 3) throw new IllegalArgumentException("Quaternary values must be 0-3");
        return Math.max(a, b);
    }

    /**
     * Quaternary Inverter (3 - x).
     * @param a Value (0-3).
     * @return 3 - a.
     */
    public static int quatInvert(int a) {
        if (a < 0 || a > 3) throw new IllegalArgumentException("Quaternary value must be 0-3");
        return 3 - a;
    }

    /**
     * Quaternary Successor (cycles: 0->1,1->2,2->3,3->0).
     * @param a Value (0-3).
     * @return (a + 1) % 4.
     */
    public static int quatSuccessor(int a) {
        if (a < 0 || a > 3) throw new IllegalArgumentException("Quaternary value must be 0-3");
        return (a + 1) % 4;
    }

    // DNA sequence operations (element-wise, safe math)
    /**
     * Applies quaternary MIN gate element-wise on two DNA sequences.
     * Sequences must be same length.
     * @param dna1 First DNA sequence.
     * @param dna2 Second DNA sequence.
     * @return New DNA sequence with MIN applied.
     */
    public static String dnaQuatMin(String dna1, String dna2) {
        validateDNA(dna1);
        validateDNA(dna2);
        if (dna1.length() != dna2.length()) {
            throw new IllegalArgumentException("DNA sequences must be same length");
        }
        StringBuilder result = new StringBuilder();
        for (int i = 0; i < dna1.length(); i++) {
            int a = DNA_TO_QUAT.get(Character.toUpperCase(dna1.charAt(i)));
            int b = DNA_TO_QUAT.get(Character.toUpperCase(dna2.charAt(i)));
            result.append(QUAT_TO_DNA.get(quatMin(a, b)));
        }
        return result.toString();
    }

    /**
     * Applies quaternary MAX gate element-wise on two DNA sequences.
     * @param dna1 First DNA sequence.
     * @param dna2 Second DNA sequence.
     * @return New DNA sequence with MAX applied.
     */
    public static String dnaQuatMax(String dna1, String dna2) {
        validateDNA(dna1);
        validateDNA(dna2);
        if (dna1.length() != dna2.length()) {
            throw new IllegalArgumentException("DNA sequences must be same length");
        }
        StringBuilder result = new StringBuilder();
        for (int i = 0; i < dna1.length(); i++) {
            int a = DNA_TO_QUAT.get(Character.toUpperCase(dna1.charAt(i)));
            int b = DNA_TO_QUAT.get(Character.toUpperCase(dna2.charAt(i)));
            result.append(QUAT_TO_DNA.get(quatMax(a, b)));
        }
        return result.toString();
    }

    /**
     * Applies quaternary inverter to each base in a DNA sequence.
     * @param dna DNA sequence.
     * @return Inverted DNA sequence.
     */
    public static String dnaQuatInvert(String dna) {
        validateDNA(dna);
        StringBuilder result = new StringBuilder();
        for (char c : dna.toCharArray()) {
            int val = DNA_TO_QUAT.get(Character.toUpperCase(c));
            result.append(QUAT_TO_DNA.get(quatInvert(val)));
        }
        return result.toString();
    }

    /**
     * Adds two DNA sequences as base-4 numbers, using BigInteger for safe arbitrary size.
     * Sequences are treated as big-endian (leftmost base is MSB).
     * @param dna1 First DNA sequence (base-4 number).
     * @param dna2 Second DNA sequence (base-4 number).
     * @return Resulting DNA sequence (sum in base-4).
     */
    public static String dnaQuatAdd(String dna1, String dna2) {
        validateDNA(dna1);
        validateDNA(dna2);
        BigInteger num1 = dnaToBigInteger(dna1);
        BigInteger num2 = dnaToBigInteger(dna2);
        BigInteger sum = num1.add(num2);
        return bigIntegerToDna(sum);
    }

    // Helper: Convert DNA sequence to BigInteger (base-4)
    private static BigInteger dnaToBigInteger(String dna) {
        BigInteger num = BigInteger.ZERO;
        BigInteger base4 = BigInteger.valueOf(4);
        for (char c : dna.toUpperCase().toCharArray()) {
            int digit = DNA_TO_QUAT.get(c);
            num = num.multiply(base4).add(BigInteger.valueOf(digit));
        }
        return num;
    }

    // Helper: Convert BigInteger to DNA sequence (base-4)
    private static String bigIntegerToDna(BigInteger num) {
        if (num.signum() < 0) {
            throw new IllegalArgumentException("Negative numbers not supported in DNA representation");
        }
        if (num.equals(BigInteger.ZERO)) {
            return "A";  // Represent 0 as 'A'
        }
        StringBuilder dna = new StringBuilder();
        BigInteger base4 = BigInteger.valueOf(4);
        while (num.compareTo(BigInteger.ZERO) > 0) {
            int remainder = num.mod(base4).intValue();
            dna.append(QUAT_TO_DNA.get(remainder));
            num = num.divide(base4);
        }
        return dna.reverse().toString();
    }

    // Base-4096 encoding using DNA (6 bases = 1 symbol, since 4^6 = 4096)
    /**
     * Encodes a DNA sequence into a base-4096 string using the fixed alphabet.
     * Pads with 'A' if length not multiple of 6.
     * @param dna DNA sequence.
     * @return Base-4096 encoded string.
     */
    public static String encodeDnaToBase4096(String dna) {
        validateDNA(dna);
        String paddedDna = dna.toUpperCase();
        int padLength = (6 - (paddedDna.length() % 6)) % 6;
        paddedDna += "A".repeat(padLength);  // Pad with A (0)
        
        StringBuilder encoded = new StringBuilder();
        for (int i = 0; i < paddedDna.length(); i += 6) {
            String group = paddedDna.substring(i, i + 6);
            BigInteger value = dnaToBigInteger(group);
            int index = value.intValueExact();  // Safe since 4^6 = 4096
            encoded.append(BASE4096_ALPHABET.get(index));
        }
        return encoded.toString();
    }

    /**
     * Decodes a base-4096 string back to DNA sequence.
     * @param encoded Base-4096 string.
     * @return Original DNA sequence (without padding).
     */
    public static String decodeBase4096ToDna(String encoded) {
        StringBuilder dna = new StringBuilder();
        for (char symbol : encoded.toCharArray()) {
            int index = BASE4096_ALPHABET.indexOf(symbol);
            if (index == -1) {
                throw new IllegalArgumentException("Invalid base-4096 symbol: " + symbol);
            }
            BigInteger value = BigInteger.valueOf(index);
            String groupDna = bigIntegerToDna(value);
            // Pad to 6 bases with leading A (0)
            while (groupDna.length() < 6) {
                groupDna = "A" + groupDna;
            }
            dna.append(groupDna);
        }
        // Remove trailing padding A's
        return dna.toString().replaceAll("A+$", "");
    }

    // Convert Float4096 to DNA representation (integer part as base-4 DNA)
    /**
     * Converts the integer part of this Float4096 to a DNA sequence (base-4 representation).
     * @return DNA sequence.
     */
    public String toDnaSequence() {
        BigInteger intPart = this.value.toBigInteger();
        return bigIntegerToDna(intPart);
    }

    /**
     * Creates a Float4096 from a DNA sequence interpreted as base-4 number.
     * @param dna DNA sequence.
     * @return Float4096 representation.
     */
    public static Float4096 fromDnaSequence(String dna) {
        BigInteger num = dnaToBigInteger(dna);
        return new Float4096(new BigDecimal(num, PRECISION_4096));
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