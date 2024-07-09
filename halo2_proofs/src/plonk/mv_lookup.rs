use crate::helpers::SerdePrimeField;

use super::{
    circuit::Expression, read_expressions_2d_vec, read_expressions_vec, write_expressions_2d_slice,
    write_expressions_slice,
};
use ff::{Field, FromUniformBytes};
use std::{
    fmt::{self, Debug},
    io,
};

pub(crate) mod prover;
pub(crate) mod verifier;

/// Degree of lookup without inputs
pub fn base_degree(table_degree: usize) -> usize {
    // let lhs_degree = table_degree + inputs_expressions_degree + 1
    // let degree = lhs_degree + 1
    std::cmp::max(3, table_degree + 2)
}

pub fn degree_with_input(base_degree: usize, input_expression_degree: usize) -> usize {
    base_degree + input_expression_degree
}

#[derive(Clone)]
pub struct Argument<F: Field> {
    pub name: &'static str,
    pub(crate) table_expressions: Vec<Expression<F>>,
    pub(crate) inputs_expressions: Vec<Vec<Expression<F>>>,
}

impl<F: Field> Debug for Argument<F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Argument")
            .field("table_expressions", &self.table_expressions)
            .field("inputs_expressions", &self.inputs_expressions)
            .finish()
    }
}

impl<F: Field> Argument<F> {
    /// Constructs a new lookup argument.
    pub fn new(name: &'static str, table: &[Expression<F>], input: &[Vec<Expression<F>>]) -> Self {
        Self {
            name,
            table_expressions: table.to_owned(),
            inputs_expressions: input.to_owned(),
        }
    }

    pub(crate) fn required_degree(&self) -> usize {
        assert!(self
            .inputs_expressions
            .iter()
            .all(|input| input.len() == self.table_expressions.len()));

        let expr_degree = |input_expressions: &Vec<Expression<F>>| {
            let mut input_degree = 0;
            for expr in input_expressions.iter() {
                input_degree = std::cmp::max(input_degree, expr.degree());
            }

            input_degree
        };

        let inputs_expressions_degree: usize =
            self.inputs_expressions.iter().map(expr_degree).sum();

        let table_degree = expr_degree(&self.table_expressions);

        /*
            φ_i(X) = f_i(X) + α
            τ(X) = t(X) + α
            LHS = τ(X) * Π(φ_i(X)) * (ϕ(gX) - ϕ(X))
                = table_degree + sum(input_degree) + 1
            RHS = τ(X) * Π(φ_i(X)) * (∑ 1/(φ_i(X)) - m(X) / τ(X))))

            deg(q(X)) = (1 - (q_last + q_blind)) * (LHS - RHS)
                 = 1 + LHS
        */

        let lhs_degree = table_degree + inputs_expressions_degree + 1;
        let degree = lhs_degree + 1;

        // 3 = phi + q_blind + table (where table is = 1)
        // + 1 for each of inputs expressions
        std::cmp::max(3 + self.inputs_expressions.len(), degree)
    }

    /// Returns input of this argument
    pub fn input_expressions(&self) -> &Vec<Vec<Expression<F>>> {
        &self.inputs_expressions
    }

    /// Returns table of this argument
    pub fn table_expressions(&self) -> &Vec<Expression<F>> {
        &self.table_expressions
    }
}

impl<F: FromUniformBytes<64>> Argument<F> {
    /// Gets the total number of bytes in the serialization of `self`
    pub(crate) fn bytes_length(&self) -> usize {
        8 + self.inputs_expressions.iter().fold(4, |acc, e_vec| {
            acc + e_vec.iter().fold(0, |acc, e| acc + e.bytes_length())
        }) + self
            .table_expressions
            .iter()
            .fold(0, |acc, e| acc + e.bytes_length())
    }
}

impl<F: SerdePrimeField + FromUniformBytes<64>> Argument<F> {
    /// Writes an argument to a buffer.
    pub fn write<W: io::Write>(&self, writer: &mut W) -> io::Result<()> {
        // NOTE(chokobole): `self.name` is not important in the sense of creating proof.
        write_expressions_2d_slice(self.inputs_expressions.as_slice(), writer)?;
        write_expressions_slice(self.table_expressions.as_slice(), writer)?;
        Ok(())
    }

    /// Reads an argument from a buffer.
    pub fn read<R: io::Read>(reader: &mut R) -> io::Result<Self> {
        Ok(Self {
            name: "",
            inputs_expressions: read_expressions_2d_vec(reader)?,
            table_expressions: read_expressions_vec(reader)?,
        })
    }
}
