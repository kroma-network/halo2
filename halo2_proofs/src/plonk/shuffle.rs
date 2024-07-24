use crate::helpers::SerdePrimeField;

use super::{circuit::Expression, read_expressions_vec, write_expressions_slice};
use ff::{Field, FromUniformBytes};
use std::{
    fmt::{self, Debug},
    io,
};

pub(crate) mod prover;
pub(crate) mod verifier;

#[derive(Clone, PartialEq, Eq)]
pub struct Argument<F: Field> {
    pub(crate) name: String,
    pub(crate) input_expressions: Vec<Expression<F>>,
    pub(crate) shuffle_expressions: Vec<Expression<F>>,
}

impl<F: Field> Debug for Argument<F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Argument")
            .field("input_expressions", &self.input_expressions)
            .field("shuffle_expressions", &self.shuffle_expressions)
            .finish()
    }
}

impl<F: Field> Argument<F> {
    /// Constructs a new shuffle argument.
    ///
    /// `shuffle` is a sequence of `(input, shuffle)` tuples.
    pub fn new<S: AsRef<str>>(name: S, shuffle: Vec<(Expression<F>, Expression<F>)>) -> Self {
        let (input_expressions, shuffle_expressions) = shuffle.into_iter().unzip();
        Argument {
            name: name.as_ref().to_string(),
            input_expressions,
            shuffle_expressions,
        }
    }

    pub(crate) fn required_degree(&self) -> usize {
        assert_eq!(self.input_expressions.len(), self.shuffle_expressions.len());

        let mut input_degree = 1;
        for expr in self.input_expressions.iter() {
            input_degree = std::cmp::max(input_degree, expr.degree());
        }
        let mut shuffle_degree = 1;
        for expr in self.shuffle_expressions.iter() {
            shuffle_degree = std::cmp::max(shuffle_degree, expr.degree());
        }

        // (1 - (l_last + l_blind)) (z(\omega X) (s(X) + \gamma) - z(X) (a(X) + \gamma))
        std::cmp::max(2 + shuffle_degree, 2 + input_degree)
    }

    /// Returns input of this argument
    pub fn input_expressions(&self) -> &Vec<Expression<F>> {
        &self.input_expressions
    }

    /// Returns table of this argument
    pub fn shuffle_expressions(&self) -> &Vec<Expression<F>> {
        &self.shuffle_expressions
    }

    /// Returns name of this argument
    pub fn name(&self) -> &str {
        &self.name
    }
}

impl<F: FromUniformBytes<64>> Argument<F> {
    /// Gets the total number of bytes in the serialization of `self`
    pub(crate) fn bytes_length(&self) -> usize {
        8 + self
            .input_expressions
            .iter()
            .fold(0, |acc, e| acc + e.bytes_length())
            + self
                .shuffle_expressions
                .iter()
                .fold(0, |acc, e| acc + e.bytes_length())
    }
}

impl<F: SerdePrimeField + FromUniformBytes<64>> Argument<F> {
    /// Writes an argument to a buffer.
    pub fn write<W: io::Write>(&self, writer: &mut W) -> io::Result<()> {
        // NOTE(chokobole): `self.name` is not important in the sense of creating proof.
        write_expressions_slice(self.input_expressions.as_slice(), writer)?;
        write_expressions_slice(self.shuffle_expressions.as_slice(), writer)?;
        Ok(())
    }

    /// Reads an argument from a buffer.
    pub fn read<R: io::Read>(reader: &mut R) -> io::Result<Self> {
        Ok(Self {
            name: "".to_string(),
            input_expressions: read_expressions_vec(reader)?,
            shuffle_expressions: read_expressions_vec(reader)?,
        })
    }
}
