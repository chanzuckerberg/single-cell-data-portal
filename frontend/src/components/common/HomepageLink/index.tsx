import Link from "next/link";
import { FC } from "react";
import { ROUTES } from "src/common/constants/routes";
import { FEATURES } from "src/common/featureFlags/features";
import { useFeatureFlag } from "src/common/hooks/useFeatureFlag";
import { Logo } from "../Logo";

export const HomepageLink: FC = () => {
  const homeRoute = useFeatureFlag(FEATURES.FILTER)
    ? ROUTES.DATASETS
    : ROUTES.HOMEPAGE;
  return (
    <Link href={homeRoute} passHref>
      <a href="passHref">
        <Logo />
      </a>
    </Link>
  );
};
