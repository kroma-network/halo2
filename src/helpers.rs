use std::io;

use crate::arithmetic::CurveAffine;

pub(crate) trait CurveRead: CurveAffine {}

impl<C: CurveAffine> CurveRead for C {}
