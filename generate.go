// Package divideandpair hosts repository-level code generation directives.
//
// The generator lives in its own module under internal/generator so that its
// dependencies (github.com/mmcloughlin/addchain) do not leak into this
// module or the anonymized artifact.
package divideandpair

//go:generate go run -C internal/generator .
