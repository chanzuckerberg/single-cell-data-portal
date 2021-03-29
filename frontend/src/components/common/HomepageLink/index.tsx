import React, { FC } from "react";
import { ROUTES } from "src/common/constants/routes";
import { Logo } from "../Logo";
import { Link } from "./style";

export const HomepageLink: FC = () => {
  return (
    <Link href={ROUTES.HOMEPAGE}>
      {/* TODO(thuang): Remove the disable once the issue resolves: https://github.com/vercel/next.js/discussions/8207 */}
      {/* eslint-disable-next-line jsx-a11y/anchor-is-valid */}
      <a>
        <Logo />
      </a>
    </Link>
  );
};
