import { AnchorButton } from "@blueprintjs/core";
import Image from "next/image";
import React from "react";
import { ActionButton as StyledActionButton } from "src/components/common/Grid/components/ActionButton/style";

interface Props {
  imageProps: StaticImageData;
  isDisabled?: boolean;
  url?: string;
}

export default function ActionButton({
  imageProps,
  isDisabled = false,
  url,
  ...props /* Spread props to allow for data-test-id and anchor button specific attributes e.g. "target". */
}: Props): JSX.Element {
  const { src, ...restImageProps } = imageProps;
  const actionIcon = <Image alt="Download" src={src} {...restImageProps} />;
  return url ? (
    <StyledActionButton
      as={AnchorButton}
      href={url}
      icon={actionIcon}
      minimal
      {...props}
      disabled={isDisabled} /* overrides props.disabled */
    />
  ) : (
    <StyledActionButton icon={actionIcon} minimal {...props} />
  );
}
