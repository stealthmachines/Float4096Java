package com.xai.float4096;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;
import java.math.RoundingMode;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Deque;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
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
 * Includes quaternary (base-4) DNA logic: 
 * - Quaternary logic gates (MIN, MAX, Inverter, Successor) for multi-valued logic.
 * - DNA sequences as quaternary representations (A=0, C=1, G=2, T=3).
 * - Safe math operations on DNA sequences (element-wise gates, addition in base-4).
 * - Base-4096 encoding/decoding using DNA sequences, where 6 DNA bases encode one base-4096 symbol (since 4^6 = 4096).
 * - Uses a fixed alphabet of 4096 unique characters for base-4096 symbols.
 * 
 * To ensure full Turing completeness, includes a DNA-based Turing Machine simulator.
 * The Turing Machine uses DNA sequences for the tape (quaternary symbols), quaternary logic for transitions,
 * and safe BigDecimal/BigInteger operations where needed to avoid overflows.
 * Supports all standard Turing operations: read/write, move left/right, state transitions, halting.
 * Can simulate any Turing-computable function, demonstrating the library's Turing completeness elegantly.
 * 
 * Creates its own native language: QuatDnaLang, a quaternary Brainfuck variant optimized for DNA and the library.
 * Programs are DNA strings, parsed in pairs for commands, executed on a quaternary tape using the library's logic gates.
 * This language is Turing-complete, self-contained, and leverages the library's quaternary DNA features for elegance.
 * Supports recursive evaluation via 'TT' command, allowing programs to execute sub-programs from the tape, vertically integrating the language with itself.
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

    // Constructors
    public Float4096(String value) {
        this.value = new BigDecimal(value, PRECISION_4096);
    }

    public Float4096(BigDecimal value) {
        this.value = Objects.requireNonNull(value).round(PRECISION_4096);
    }

    public Float4096(double value) {
        this.value = new BigDecimal(value, PRECISION_4096);
    }

    // Arithmetic operations
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

    // Comparison and utilities
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

    // Binary logic gates using subtraction
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

    // Deterministic generation of 4096 unique characters
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

    // Quaternary DNA Logic Section

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
    public static int quatMin(int a, int b) {
        if (a < 0 || a > 3 || b < 0 || b > 3) throw new IllegalArgumentException("Quaternary values must be 0-3");
        return Math.min(a, b);
    }

    public static int quatMax(int a, int b) {
        if (a < 0 || a > 3 || b < 0 || b > 3) throw new IllegalArgumentException("Quaternary values must be 0-3");
        return Math.max(a, b);
    }

    public static int quatInvert(int a) {
        if (a < 0 || a > 3) throw new IllegalArgumentException("Quaternary value must be 0-3");
        return 3 - a;
    }

    public static int quatSuccessor(int a) {
        if (a < 0 || a > 3) throw new IllegalArgumentException("Quaternary value must be 0-3");
        return (a + 1) % 4;
    }

    public static int quatPredecessor(int a) {
        if (a < 0 || a > 3) throw new IllegalArgumentException("Quaternary value must be 0-3");
        return (a + 3) % 4;
    }

    // DNA sequence operations (element-wise, safe math)
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

    public static String dnaQuatInvert(String dna) {
        validateDNA(dna);
        StringBuilder result = new StringBuilder();
        for (char c : dna.toCharArray()) {
            int val = DNA_TO_QUAT.get(Character.toUpperCase(c));
            result.append(QUAT_TO_DNA.get(quatInvert(val)));
        }
        return result.toString();
    }

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
    public String toDnaSequence() {
        BigInteger intPart = this.value.toBigInteger();
        return bigIntegerToDna(intPart);
    }

    public static Float4096 fromDnaSequence(String dna) {
        BigInteger num = dnaToBigInteger(dna);
        return new Float4096(new BigDecimal(num, PRECISION_4096));
    }

    // DNA-based Turing Machine Simulator for full Turing completeness
    public static final class DnaTuringMachine {
        private final Map<String, Map<Character, Transition>> transitions;
        private final char blankSymbol = 'A';
        private List<Character> tape;
        private int headPosition;
        private String currentState;
        private final String haltState;

        // Transition record
        public record Transition(String newState, char newSymbol, char direction) {
            public Transition {
                if (direction != 'L' && direction != 'R') {
                    throw new IllegalArgumentException("Direction must be 'L' or 'R'");
                }
            }
        }

        private DnaTuringMachine(Builder builder) {
            this.transitions = new HashMap<>(builder.transitions);
            this.tape = new ArrayList<>();
            for (char c : builder.initialTape.toUpperCase().toCharArray()) {
                this.tape.add(c);
            }
            this.headPosition = builder.initialHeadPosition;
            this.currentState = builder.initialState;
            this.haltState = builder.haltState;
        }

        public String run(int maxSteps) {
            int step = 0;
            while (!currentState.equals(haltState) && step < maxSteps) {
                if (headPosition < 0) {
                    tape.add(0, blankSymbol);
                    headPosition = 0;
                } else if (headPosition >= tape.size()) {
                    tape.add(blankSymbol);
                }

                char currentSymbol = tape.get(headPosition);
                Map<Character, Transition> stateTrans = transitions.get(currentState);
                if (stateTrans == null) {
                    throw new RuntimeException("No transition for state: " + currentState);
                }
                Transition trans = stateTrans.get(currentSymbol);
                if (trans == null) {
                    throw new RuntimeException("No transition for state " + currentState + " and symbol " + currentSymbol);
                }

                tape.set(headPosition, trans.newSymbol);

                if (trans.direction == 'L') {
                    headPosition--;
                } else {
                    headPosition++;
                }

                currentState = trans.newState;

                step++;
            }

            if (step >= maxSteps) {
                throw new RuntimeException("Max steps exceeded; possible infinite loop");
            }

            StringBuilder finalTape = new StringBuilder();
            boolean started = false;
            for (char c : tape) {
                if (c != blankSymbol || started) {
                    started = true;
                    finalTape.append(c);
                }
            }
            String tapeStr = finalTape.toString().replaceAll(blankSymbol + "+$", "");
            return tapeStr.isEmpty() ? String.valueOf(blankSymbol) : tapeStr;
        }

        public static final class Builder {
            private final Map<String, Map<Character, Transition>> transitions = new HashMap<>();
            private String initialTape = "";
            private int initialHeadPosition = 0;
            private String initialState = "start";
            private String haltState = "halt";

            public Builder addTransition(String state, char symbol, String newState, char newSymbol, char direction) {
                transitions.computeIfAbsent(state, k -> new HashMap<>()).put(Character.toUpperCase(symbol), new Transition(newState, Character.toUpperCase(newSymbol), direction));
                return this;
            }

            public Builder initialTape(String tape) {
                validateDNA(tape);
                this.initialTape = tape;
                return this;
            }

            public Builder initialHeadPosition(int position) {
                if (position < 0) throw new IllegalArgumentException("Head position must be non-negative");
                this.initialHeadPosition = position;
                return this;
            }

            public Builder initialState(String state) {
                this.initialState = state;
                return this;
            }

            public Builder haltState(String state) {
                this.haltState = state;
                return this;
            }

            public DnaTuringMachine build() {
                if (transitions.isEmpty()) {
                    throw new IllegalStateException("No transitions defined");
                }
                return new DnaTuringMachine(this);
            }
        }

        public static DnaTuringMachine createIncrementer(String initialDna) {
            Builder builder = new Builder().initialTape(initialDna).initialState("inc").haltState("halt").initialHeadPosition(initialDna.length() - 1);
            for (char symbol : "ACGT".toCharArray()) {
                int val = DNA_TO_QUAT.get(symbol);
                int newVal = quatSuccessor(val);
                char newSymbol = QUAT_TO_DNA.get(newVal);
                if (newVal == 0) {
                    // Carry
                    builder.addTransition("inc", symbol, "inc", 'A', 'L');
                } else {
                    // No carry
                    builder.addTransition("inc", symbol, "halt", newSymbol, 'R');
                }
            }
            // If carry at left end, extend with 'C' (1)
            builder.addTransition("inc", 'A', "halt", 'C', 'R');
            return builder.build();
        }
    }

    // Native Language: QuatDnaLang - A quaternary Brainfuck variant optimized for DNA and the library
    public static final class QuatDnaLangInterpreter {
        private static final Map<String, Character> COMMANDS = new HashMap<>();
        static {
            COMMANDS.put("AA", '+'); // Successor
            COMMANDS.put("AC", '-'); // Predecessor
            COMMANDS.put("AG", '>'); // Right
            COMMANDS.put("AT", '<'); // Left
            COMMANDS.put("CA", '.'); // Output
            COMMANDS.put("CC", ','); // Input
            COMMANDS.put("CG", '['); // Loop start
            COMMANDS.put("CT", ']'); // Loop end
            COMMANDS.put("GA", 'I'); // Invert
            COMMANDS.put("GC", 'M'); // Min with next
            COMMANDS.put("GG", 'X'); // Max with next
            COMMANDS.put("TT", 'E'); // Eval subprogram
        }

        public static String execute(String program, String input, int maxSteps) {
            validateDNA(program);
            if (program.length() % 2 != 0) {
                throw new IllegalArgumentException("Program length must be even for paired commands");
            }
            validateDNA(input);

            // Parse program to command list
            List<Character> instructions = new ArrayList<>();
            for (int i = 0; i < program.length(); i += 2) {
                String pair = program.substring(i, i + 2).toUpperCase();
                Character cmd = COMMANDS.get(pair);
                if (cmd == null) continue; // Ignore invalid as nop
                instructions.add(cmd);
            }

            // Build bracket map for loops
            Map<Integer, Integer> bracketMap = new HashMap<>();
            Deque<Integer> stack = new LinkedList<>();
            for (int i = 0; i < instructions.size(); i++) {
                char cmd = instructions.get(i);
                if (cmd == '[') {
                    stack.push(i);
                } else if (cmd == ']') {
                    if (stack.isEmpty()) {
                        throw new IllegalArgumentException("Mismatched brackets");
                    }
                    int open = stack.pop();
                    bracketMap.put(open, i);
                    bracketMap.put(i, open);
                }
            }
            if (!stack.isEmpty()) {
                throw new IllegalArgumentException("Mismatched brackets");
            }

            // Tape: quaternary cells
            List<Integer> tape = new ArrayList<>();
            tape.add(0); // Initial cell
            int head = 0;

            // Input as list of ints
            List<Integer> inputList = new ArrayList<>();
            for (char c : input.toUpperCase().toCharArray()) {
                inputList.add(DNA_TO_QUAT.get(c));
            }
            int inputIndex = 0;

            // Output
            StringBuilder output = new StringBuilder();

            // Run
            int pc = 0; // Program counter
            int step = 0;
            while (pc < instructions.size() && step < maxSteps) {
                char cmd = instructions.get(pc);
                switch (cmd) {
                    case '+' -> tape.set(head, quatSuccessor(tape.get(head)));
                    case '-' -> tape.set(head, quatPredecessor(tape.get(head)));
                    case '>' -> {
                        head++;
                        if (head >= tape.size()) tape.add(0);
                    }
                    case '<' -> {
                        head--;
                        if (head < 0) {
                            tape.add(0, 0);
                            head = 0;
                        }
                    }
                    case '.' -> output.append(QUAT_TO_DNA.get(tape.get(head)));
                    case ',' -> {
                        if (inputIndex < inputList.size()) {
                            tape.set(head, inputList.get(inputIndex++));
                        } else {
                            tape.set(head, 0);
                        }
                    }
                    case '[' -> {
                        if (tape.get(head) == 0) pc = bracketMap.get(pc);
                    }
                    case ']' -> {
                        if (tape.get(head) != 0) pc = bracketMap.get(pc);
                    }
                    case 'I' -> tape.set(head, quatInvert(tape.get(head)));
                    case 'M' -> {
                        if (head + 1 >= tape.size()) tape.add(0);
                        tape.set(head, quatMin(tape.get(head), tape.get(head + 1)));
                    }
                    case 'X' -> {
                        if (head + 1 >= tape.size()) tape.add(0);
                        tape.set(head, quatMax(tape.get(head), tape.get(head + 1)));
                    }
                    case 'E' -> {
                        head++;
                        StringBuilder subProgram = new StringBuilder();
                        while (true) {
                            if (head + 1 >= tape.size()) break;
                            int a = tape.get(head);
                            int b = tape.get(head + 1);
                            if (a == 0 && b == 0) break;
                            subProgram.append(QUAT_TO_DNA.get(a)).append(QUAT_TO_DNA.get(b));
                            head += 2;
                        }
                        String subOutput = execute(subProgram.toString(), "", maxSteps - step);
                        // Put subOutput back on tape starting from current head
                        head -= subProgram.length(); // Reset head to start of subprogram area
                        for (int i = 0; i < subOutput.length(); i++) {
                            if (head + i >= tape.size()) tape.add(0);
                            tape.set(head + i, DNA_TO_QUAT.get(subOutput.charAt(i)));
                        }
                        head += subOutput.length() - 1; // Move head to end of output
                    }
                }
                pc++;
                step++;
            }

            if (step >= maxSteps) {
                throw new RuntimeException("Max steps exceeded; possible infinite loop");
            }

            return output.toString();
        }
    }

    // Constant calculations
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