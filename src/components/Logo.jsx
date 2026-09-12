export function LogoMark({ size = 28, className = "" }) {
  return (
    <svg
      width={size}
      height={size}
      viewBox="0 0 32 32"
      fill="none"
      xmlns="http://www.w3.org/2000/svg"
      className={className}
      aria-hidden="true"
    >
      <polyline
        points="9,6 23,13 9,20 23,27"
        stroke="#E0576B"
        strokeWidth="1.8"
        strokeLinecap="square"
        strokeLinejoin="miter"
      />
      <polyline
        points="23,6 9,13 23,20 9,27"
        stroke="#F6F7F9"
        strokeWidth="1.8"
        strokeLinecap="square"
        strokeLinejoin="miter"
      />
    </svg>
  );
}

export default function Logo({ size = 26, className = "" }) {
  return (
    <div className={`flex items-center gap-2.5 ${className}`}>
      <LogoMark size={size} />
      <span className="text-[20px] font-extrabold tracking-tight">
        Apta<span className="text-accent">Scope</span>
      </span>
    </div>
  );
}
