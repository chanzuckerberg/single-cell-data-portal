import React, { FC } from "react";
import { Logo } from "../Logo";
import { Link } from "./style";

interface Props {
  small?: boolean;
}

export const HomepageLink: FC<Props> = ({ small }) => {
  return (
    <Link to="/">
      <Logo small={small} />
    </Link>
  );
};
