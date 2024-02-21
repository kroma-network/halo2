use crate::helpers::SerdePrimeField;

use super::{circuit::Expression, read_expressions_vec, write_expressions_slice};
use ff::Field;
use std::{
    fmt::{self, Debug},
    io,
};

pub(crate) mod prover;
pub(crate) mod verifier;

#[derive(Clone)]
pub struct Argument<F: Field> {
    pub name: &'static str,
    pub input_expressions: Vec<Expression<F>>,
    pub table_expressions: Vec<Expression<F>>,
}

impl<F: Field> Debug for Argument<F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Argument")
            .field("input_expressions", &self.input_expressions)
            .field("table_expressions", &self.table_expressions)
            .finish()
    }
}

impl<F: Field> Argument<F> {
    /// Constructs a new lookup argument.
    ///
    /// `table_map` is a sequence of `(input, table)` tuples.
    pub fn new(name: &'static str, table_map: Vec<(Expression<F>, Expression<F>)>) -> Self {
        let (input_expressions, table_expressions) = table_map.into_iter().unzip();
        Argument {
            name,
            input_expressions,
            table_expressions,
        }
    }

    pub(crate) fn required_degree(&self) -> usize {
        assert_eq!(self.input_expressions.len(), self.table_expressions.len());

        // The first value in the permutation poly should be one.
        // degree 2:
        // l_0(X) * (1 - z(X)) = 0
        //
        // The "last" value in the permutation poly should be a boolean, for
        // completeness and soundness.
        // degree 3:
        // l_last(X) * (z(X)^2 - z(X)) = 0
        //
        // Enable the permutation argument for only the rows involved.
        // degree (2 + input_degree + table_degree) or 4, whichever is larger:
        // (1 - (l_last(X) + l_blind(X))) * (
        //   z(\omega X) (a'(X) + \beta) (s'(X) + \gamma)
        //   - z(X) (\theta^{m-1} a_0(X) + ... + a_{m-1}(X) + \beta) (\theta^{m-1} s_0(X) + ... + s_{m-1}(X) + \gamma)
        // ) = 0
        //
        // The first two values of a' and s' should be the same.
        // degree 2:
        // l_0(X) * (a'(X) - s'(X)) = 0
        //
        // Either the two values are the same, or the previous
        // value of a' is the same as the current value.
        // degree 3:
        // (1 - (l_last(X) + l_blind(X))) * (a′(X) − s′(X))⋅(a′(X) − a′(\omega^{-1} X)) = 0
        let mut input_degree = 1;
        for expr in self.input_expressions.iter() {
            input_degree = std::cmp::max(input_degree, expr.degree());
        }
        let mut table_degree = 1;
        for expr in self.table_expressions.iter() {
            table_degree = std::cmp::max(table_degree, expr.degree());
        }

        // In practice because input_degree and table_degree are initialized to
        // one, the latter half of this max() invocation is at least 4 always,
        // rendering this call pointless except to be explicit in case we change
        // the initialization of input_degree/table_degree in the future.
        std::cmp::max(
            // (1 - (l_last + l_blind)) z(\omega X) (a'(X) + \beta) (s'(X) + \gamma)
            4,
            // (1 - (l_last + l_blind)) z(X) (\theta^{m-1} a_0(X) + ... + a_{m-1}(X) + \beta) (\theta^{m-1} s_0(X) + ... + s_{m-1}(X) + \gamma)
            2 + input_degree + table_degree,
        )
    }

    /// Returns input of this argument
    pub fn input_expressions(&self) -> &Vec<Expression<F>> {
        &self.input_expressions
    }

    /// Returns table of this argument
    pub fn table_expressions(&self) -> &Vec<Expression<F>> {
        &self.table_expressions
    }
}

impl<F: SerdePrimeField> Argument<F> {
    /// Gets the total number of bytes in the serialization of `self`
    pub(crate) fn bytes_length(&self) -> usize {
        8 + self
            .input_expressions
            .iter()
            .fold(0, |acc, e| acc + e.bytes_length())
            + self
                .table_expressions
                .iter()
                .fold(0, |acc, e| acc + e.bytes_length())
    }

    /// Writes an argument to a buffer.
    pub fn write<W: io::Write>(&self, writer: &mut W) -> io::Result<()> {
        // NOTE(chokobole): `self.name` is not important in the sense of creating proof.
        write_expressions_slice(self.input_expressions.as_slice(), writer)?;
        write_expressions_slice(self.table_expressions.as_slice(), writer)?;
        Ok(())
    }

    /// Reads an argument from a buffer.
    pub fn read<R: io::Read>(reader: &mut R) -> io::Result<Self> {
        Ok(Self {
            name: "",
            input_expressions: read_expressions_vec(reader)?,
            table_expressions: read_expressions_vec(reader)?,
        })
    }
}
