import Link from "next/link";
import React from "react";
import { StyledAnchor } from "src/components/common/Grid/components/LinkCell/style";

interface Props {
  url: string;
  value: string;
}

export default function LinkCell({
  url,
  value,
  ...props /* Spread props to allow for data-testid. */
}: Props): JSX.Element {
  return (
    <Link href={url} passHref>
      <StyledAnchor href="passHref" {...props}>
        {value}
      </StyledAnchor>
    </Link>
  );
}
